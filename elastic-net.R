# =============================================================================
# Elastic Net Regression: Predicting MPRA Bias Scores from TF Motif Features
#
# This script trains elastic net regression models (via glmnet) to predict
# two types of fragment-level bias scores derived from a massively parallel
# reporter assay (MPRA):
#
#   1. Promoter bias  – z-score reflecting higher activity with INS vs. SCP1
#   2. Position bias  – z-score reflecting higher activity upstream vs. downstream
#
# For each outcome, the script:
#   - Loads pre-computed fragment-level z-scores and TF motif overlap scores
#   - Subsets features to a curated motif list
#   - Performs an 80/20 train/test split
#   - Fits cross-validated elastic net models across a grid of alpha values
#     (alpha = 0 is ridge, alpha = 1 is LASSO, values between are elastic net)
#   - Selects the best alpha by test-set MSE
#   - Reports model performance (R², MSE, RMSE, Pearson r)
#   - Saves diagnostic observed-vs-predicted scatter plots and a metrics table
#
# Input files (relative to working directory):
#   data/motif_list.encode_1KG.trimmed.txt   – curated motif names (one per line)
#   data/wald_prom.zscores.tsv               – promoter bias z-scores + motif features
#   data/wald_pos.zscores.tsv                – position bias z-scores + motif features
#
# Outputs:
#   output/metrics.tsv                       – train/test performance metrics
#   output/prom_test_plot.png                – promoter model obs-vs-pred (full range)
#   output/prom_test_plot_zoom.png           – promoter model obs-vs-pred (zoomed)
#   output/pos_test_plot.png                 – position model obs-vs-pred (full range)
#   output/pos_test_plot_zoom.png            – position model obs-vs-pred (zoomed)
# =============================================================================

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(glmnet))

set.seed(2022)

# -----------------------------------------------------------------------------
# Helper functions
# -----------------------------------------------------------------------------

#' Compute regression performance metrics and print a summary line.
#'
#' @param y_obs   Numeric vector of observed values.
#' @param y_pred  Numeric vector of predicted values (same length as y_obs).
#' @param label   Character string used as a prefix in the printed output.
#' @return Invisible list with named elements: r2, mse, rmse, pearson_r, pearson_p.
compute_metrics <- function(y_obs, y_pred, label = "") {
  residuals <- y_obs - y_pred
  ss_res    <- sum(residuals^2)
  ss_tot    <- sum((y_obs - mean(y_obs))^2)
  r2        <- 1 - ss_res / ss_tot
  mse       <- mean(residuals^2)
  rmse      <- sqrt(mse)

  pearson   <- cor.test(y_obs, y_pred, method = "pearson")
  pearson_r <- pearson$estimate
  pearson_p <- pearson$p.value

  cat(sprintf(
    "[%s]  R² = %.4f  |  MSE = %.4f  |  RMSE = %.4f  |  Pearson r = %.4f  |  Pearson p = %.4e\n",
    label, r2, mse, rmse, pearson_r, pearson_p
  ))

  invisible(list(
    r2        = r2,
    mse       = mse,
    rmse      = rmse,
    pearson_r = pearson_r,
    pearson_p = pearson_p
  ))
}


#' Build an observed-vs-predicted scatter plot for a regression model.
#'
#' A 1:1 identity line (red dashed) and zero-intercept cross-hairs are drawn
#' for reference. When annotate_plot = TRUE, R² and Pearson r (with significance
#' indicator) are added as a text annotation in the upper-left.
#'
#' @param y_obs          Numeric vector of observed z-scores.
#' @param y_pred         Numeric vector of predicted z-scores.
#' @param metrics        List returned by compute_metrics().
#' @param title_str      Plot title string.
#' @param annotate_plot  Logical; whether to add the R²/r annotation.
#' @return A ggplot object.
make_obs_pred_plot <- function(y_obs, y_pred, metrics, title_str,
                               annotate_plot = TRUE) {
  df <- data.frame(obs = y_obs, pred = y_pred)

  annotation_label <- ""
  if (annotate_plot) {
    sig_str <- if (metrics$pearson_p < 0.05) "p < 0.05" else "p > 0.05"
    annotation_label <- sprintf("R² = %.3f\nr = %.3f, %s",
                                metrics$r2, metrics$pearson_r, sig_str)
  }

  x_pos <- min(df$obs)  + 0.05 * diff(range(df$obs))
  y_pos <- max(df$pred) - 0.05 * diff(range(df$pred))

  ggplot(df, aes(x = obs, y = pred)) +
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = 0) +
    geom_point(alpha = 0.3, size = 0.8, colour = "steelblue") +
    geom_abline(slope = 1, intercept = 0, colour = "red", linetype = "dashed") +
    annotate("text", x = x_pos, y = y_pos, label = annotation_label,
             hjust = 0, vjust = 1, size = 3.5) +
    labs(title = title_str,
         x     = "Observed bias score",
         y     = "Predicted bias score") +
    theme_minimal(base_size = 11) +
    theme(panel.background = element_blank(),
          panel.grid       = element_blank())
}


#' Fit cross-validated elastic net models across an alpha grid.
#'
#' Trains one cv.glmnet model per alpha value (0, 0.05, 0.10, …, 1.0) and
#' returns both the list of fitted models and a data frame of test-set MSE
#' for each alpha.
#'
#' @param x_train  Training feature matrix.
#' @param y_train  Training response vector.
#' @param x_test   Test feature matrix.
#' @param y_test   Test response vector.
#' @return List with elements:
#'   \describe{
#'     \item{models}{Named list of cv.glmnet objects, keyed by "alpha<value>".}
#'     \item{alpha_mse}{Data frame with columns alpha, mse, model_name.}
#'   }
fit_elastic_net_grid <- function(x_train, y_train, x_test, y_test) {
  alpha_values <- seq(0, 1, by = 0.05)
  models       <- list()
  alpha_mse    <- data.frame()

  for (alpha in alpha_values) {
    model_name         <- paste0("alpha", alpha)
    models[[model_name]] <- cv.glmnet(
      x_train, y_train,
      type.measure = "mse",
      alpha        = alpha,
      family       = "gaussian"
    )

    predicted <- predict(models[[model_name]],
                         s    = models[[model_name]]$lambda.min,
                         newx = x_test)
    mse       <- mean((y_test - predicted)^2)
    alpha_mse <- rbind(alpha_mse,
                       data.frame(alpha = alpha, mse = mse, model_name = model_name))
  }

  list(models = models, alpha_mse = alpha_mse)
}


#' Extract non-zero, non-intercept coefficients from a fitted cv.glmnet model.
#'
#' @param model  A cv.glmnet object.
#' @return Data frame with columns motif and coefficient.
extract_nonzero_coefs <- function(model) {
  coef_df        <- as.data.frame(as.matrix(coef(model, s = "lambda.min")))
  coef_df$motif  <- rownames(coef_df)
  rownames(coef_df) <- NULL
  colnames(coef_df)[1] <- "coefficient"
  subset(coef_df, coefficient != 0 & motif != "(Intercept)")
}

# -----------------------------------------------------------------------------
# Load inputs
# -----------------------------------------------------------------------------

# Curated TF motif list (one motif name per line); z-score column is also kept
motif_list <- c("zscore",
                scan("data/motif_list.encode_1KG.trimmed.txt", character()))

# Fragment-level z-scores and motif overlap scores
prom_df <- read_tsv("data/wald_prom.zscores.tsv")   # promoter bias
pos_df  <- read_tsv("data/wald_pos.zscores.tsv")    # position bias

# Subset both data frames to the curated motif columns
prom_features <- select(prom_df, which(names(prom_df) %in% motif_list))
pos_features  <- select(pos_df,  which(names(pos_df)  %in% motif_list))

# -----------------------------------------------------------------------------
# Train / test split (shared 80/20 split indices)
# Note: a single split is used so that promoter and position models are
# evaluated on the same held-out fragments.
# -----------------------------------------------------------------------------

n_fragments <- nrow(prom_df)
train_idx   <- sample(seq_len(n_fragments), size = 0.8 * n_fragments, replace = FALSE)

# Design matrices (model.matrix drops the intercept column via [,-1])
prom_x <- model.matrix(zscore ~ ., prom_features)[, -1]
prom_y <- prom_features$zscore

pos_x  <- model.matrix(zscore ~ ., pos_features)[, -1]
pos_y  <- pos_features$zscore

prom_x_train <- prom_x[train_idx, ];   prom_x_test <- prom_x[-train_idx, ]
prom_y_train <- prom_y[train_idx];     prom_y_test <- prom_y[-train_idx]

pos_x_train  <- pos_x[train_idx, ];    pos_x_test  <- pos_x[-train_idx, ]
pos_y_train  <- pos_y[train_idx];      pos_y_test  <- pos_y[-train_idx]

# -----------------------------------------------------------------------------
# Fit elastic net models across alpha grid
# -----------------------------------------------------------------------------

cat("Fitting promoter bias models...\n")
prom_fit <- fit_elastic_net_grid(prom_x_train, prom_y_train,
                                 prom_x_test,  prom_y_test)

cat("Fitting position bias models...\n")
pos_fit  <- fit_elastic_net_grid(pos_x_train,  pos_y_train,
                                 pos_x_test,   pos_y_test)

# Select the alpha with the lowest test-set MSE
best_prom_alpha <- prom_fit$alpha_mse[which.min(prom_fit$alpha_mse$mse), ]$alpha
best_pos_alpha  <- pos_fit$alpha_mse[which.min(pos_fit$alpha_mse$mse),  ]$alpha

cat(sprintf("Best promoter alpha: %.2f\n", best_prom_alpha))
cat(sprintf("Best position alpha: %.2f\n", best_pos_alpha))

final_prom_model <- prom_fit$models[[paste0("alpha", best_prom_alpha)]]
final_pos_model  <- pos_fit$models[[paste0("alpha",  best_pos_alpha)]]

# -----------------------------------------------------------------------------
# Model performance metrics
# -----------------------------------------------------------------------------

prom_train_pred <- as.numeric(predict(final_prom_model, s = "lambda.min", newx = prom_x_train))
prom_test_pred  <- as.numeric(predict(final_prom_model, s = "lambda.min", newx = prom_x_test))

pos_train_pred  <- as.numeric(predict(final_pos_model,  s = "lambda.min", newx = pos_x_train))
pos_test_pred   <- as.numeric(predict(final_pos_model,  s = "lambda.min", newx = pos_x_test))

prom_train_metrics <- compute_metrics(prom_y_train, prom_train_pred, "Promoter Train")
prom_test_metrics  <- compute_metrics(prom_y_test,  prom_test_pred,  "Promoter Test ")
pos_train_metrics  <- compute_metrics(pos_y_train,  pos_train_pred,  "Position Train")
pos_test_metrics   <- compute_metrics(pos_y_test,   pos_test_pred,   "Position Test ")

cat("\n===== Variance explained (test set) =====\n")
cat(sprintf("Promoter: explained = %.1f%%,  unexplained = %.1f%%\n",
            prom_test_metrics$r2 * 100, (1 - prom_test_metrics$r2) * 100))
cat(sprintf("Position: explained = %.1f%%,  unexplained = %.1f%%\n",
            pos_test_metrics$r2  * 100, (1 - pos_test_metrics$r2)  * 100))

# Save metrics table
metrics_df <- tribble(
  ~model,      ~set,    ~r2,                    ~mse,                   ~rmse,                   ~pearson_r,                  ~pearson_p,
  "Promoter",  "Train", prom_train_metrics$r2,  prom_train_metrics$mse, prom_train_metrics$rmse, prom_train_metrics$pearson_r, prom_train_metrics$pearson_p,
  "Promoter",  "Test",  prom_test_metrics$r2,   prom_test_metrics$mse,  prom_test_metrics$rmse,  prom_test_metrics$pearson_r,  prom_test_metrics$pearson_p,
  "Position",  "Train", pos_train_metrics$r2,   pos_train_metrics$mse,  pos_train_metrics$rmse,  pos_train_metrics$pearson_r,  pos_train_metrics$pearson_p,
  "Position",  "Test",  pos_test_metrics$r2,    pos_test_metrics$mse,   pos_test_metrics$rmse,   pos_test_metrics$pearson_r,   pos_test_metrics$pearson_p
)
write_tsv(metrics_df, "output/metrics.tsv")

# -----------------------------------------------------------------------------
# Non-zero model coefficients
# -----------------------------------------------------------------------------

prom_nonzero_coefs <- extract_nonzero_coefs(final_prom_model)
pos_nonzero_coefs  <- extract_nonzero_coefs(final_pos_model)

cat(sprintf("\nPromoter model: %d non-zero motif coefficients\n", nrow(prom_nonzero_coefs)))
cat(sprintf("Position model: %d non-zero motif coefficients\n",  nrow(pos_nonzero_coefs)))

# -----------------------------------------------------------------------------
# Diagnostic plots: observed vs. predicted scatter
# -----------------------------------------------------------------------------

# Full-range plots with highlighted region of interest (dashed rectangle)
p_prom_test <- make_obs_pred_plot(
    prom_y_test, prom_test_pred, prom_test_metrics, "Promoter – Test"
  ) +
  annotate("rect",
           xmin = 0, xmax = 4, ymin = 0, ymax = 4,
           fill = NA, color = "black", linetype = "dashed", linewidth = 0.8)

p_pos_test <- make_obs_pred_plot(
    pos_y_test, pos_test_pred, pos_test_metrics, "Position – Test"
  ) +
  annotate("rect",
           xmin = 0, xmax = 40, ymin = 0, ymax = 40,
           fill = NA, color = "black", linetype = "dashed", linewidth = 0.8)

# Zoomed plots focusing on the highlighted region
p_prom_test_zoom <- make_obs_pred_plot(
    prom_y_test, prom_test_pred, prom_test_metrics, "", annotate_plot = FALSE
  ) + xlim(0, 4.5) + ylim(0, 4.5)

p_pos_test_zoom <- make_obs_pred_plot(
    pos_y_test, pos_test_pred, pos_test_metrics, "", annotate_plot = FALSE
  ) + xlim(0, 40) + ylim(0, 40)

# Save plots
ggsave("output/prom_test_plot.png",      plot = p_prom_test,      width = 2.5, height = 2.5, dpi = 300)
ggsave("output/pos_test_plot.png",       plot = p_pos_test,       width = 2.5, height = 2.5, dpi = 300)
ggsave("output/prom_test_plot_zoom.png", plot = p_prom_test_zoom, width = 2.5, height = 2.5, dpi = 300)
ggsave("output/pos_test_plot_zoom.png",  plot = p_pos_test_zoom,  width = 2.5, height = 2.5, dpi = 300)