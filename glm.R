# =============================================================================
# MPRA Fragment Activity GLM Analysis
# =============================================================================
# Loads per-fragment activity z-scores from MPRA
# annotates fragments with ChromHMM states (islet promoters/enhancers),
# ATAC-seq peak overlap, and CAGE-strand orientation. Produces:
#   1. Bar charts of annotation/strand counts across fragments
#   2. Raincloud plots of z-score distributions by strand concordance & element type
#   3. Univariate GLMs predicting fragment activity from each annotation
#   4. Model diagnostic plots (residuals vs fitted, Q-Q, scale-location, histogram)
#   5. Per-configuration GLMs for each of the 4 plasmid configurations tested
#
# Expected directory layout (relative to this script):
#   data/
#     quant.tsv                          - MPRA per-fragment z-scores
#     oligo_coords_with_strand.tsv       - Fragment coordinates + strand + sequence
#     cage_tissue_overlap/               - ChromHMM/ATAC overlap files (*_loj198bp.tsv)
#     data/config_activities/            - Per-configuration activity files
#   figures/                             - Output directory for saved plots
# =============================================================================
 
suppressPackageStartupMessages({
  library(tidyverse)
  library(glmnet)
  library(cowplot)
  library(glue)
  library(ggdist)   
})

set.seed(2023)
 
# =============================================================================
# 1. Load MPRA z-scores and parse fragment coordinates
# =============================================================================

# Load per-fragment activity z-scores. Row names encode fragment coordinates
# as "TC_chr_start_end".
fragment_zscores <- read.table("data/quant.tsv", header = TRUE)

# Move row names into a column, then split into coordinate fields
fragment_zscores$coord <- rownames(fragment_zscores)
rownames(fragment_zscores) <- NULL

fragment_zscores <- separate_wider_delim(
  fragment_zscores,
  cols  = coord,
  delim = "_",
  names = c("TC", "chr", "start", "end")
)

fragment_coords <- fragment_zscores %>%
  select(zscore, chr, start, end) %>%
  mutate(
    start = as.numeric(start),
    end   = as.numeric(end)
  )
 
 
# =============================================================================
# 2. Add CAGE strand orientation
# =============================================================================

# Load fragment coordinates with strand assignment derived from CAGE data.
# Columns: chr, start, end, chr2, start2, end2, sequence, strand
strand_annotations <- suppressMessages(
  read_tsv(
    "data/oligo_coords_with_strand.tsv",
    col_names = c("chr", "start", "end", "chr2", "start2", "end2", "sequence", "strand")
  )
) %>%
  select(chr, start, end, strand, sequence)
 
# =============================================================================
# 3. Add ChromHMM and ATAC-seq overlap annotations
# =============================================================================
 
# Each overlap file encodes the number of bp of a fragment overlapping a given
# ChromHMM state or ATAC-seq peak set (suffix: *_loj198bp.tsv).
# File naming convention: <prefix>_<prefix2>_<tissue>_<annotation>_loj198bp.tsv
# We restrict to "islet" tissue files.
tissue_overlap_files <- list.files(
  path     = "data/cage_tissue_overlap",
  pattern  = "*_loj198bp.tsv",
  full.names = TRUE
)
tissue_overlap_files <- tissue_overlap_files[grepl("islet", tissue_overlap_files)]
 
# For each overlap file, create a binary annotation column (0/1 overlap)
# and attach it to the fragment coordinate data frame.
annotation_cols <- lapply(tissue_overlap_files, function(file) {
 
  # Parse tissue and annotation type from filename
  file_parts <- str_split(basename(file), "_")[[1]]
  tissue      <- file_parts[3]
  annot_type  <- file_parts[4]
 
  # Load overlap counts (chr, start, end, overlap_count)
  overlap_counts <- suppressMessages(
    read_tsv(file, col_names = c("chr", "start", "end", "overlap_count"))
  ) %>%
    select(chr, start, end, overlap_count)
 
  # Join to fragment coordinates and recode as binary overlap indicator
  fragment_annot <- inner_join(fragment_coords, overlap_counts,
                               by = c("chr", "start", "end")) %>%
    mutate(overlap_binary = as.integer(overlap_count > 0))
 
  # Construct column name from tissue and annotation type
  col_name <- if (annot_type == "promoters") {
    tolower(paste(tissue, "promoter", sep = "_"))
  } else if (annot_type == "enhancers") {
    tolower(paste(tissue, "enhancer", sep = "_"))
  } else {
    tolower(paste(tissue, "atac", sep = "_"))
  }
 
  select(rename(fragment_annot, !!col_name := overlap_binary), all_of(col_name))
})
 
# Combine all annotation columns and join with strand data
fragment_annotations <- do.call(cbind, annotation_cols)
fragment_annotated   <- cbind(fragment_strand, fragment_annotations) %>%
  select(-sequence)   # sequence not needed for modelling
 
# Re-attach z-scores (left_join ensures z-scores align to annotated fragments)
fragment_annotated <- left_join(fragment_annotated, fragment_coords,
                                by = c("chr", "start", "end"))
 
# Add a unique fragment ID column matching the format used in per-config files
fragment_annotated <- fragment_annotated %>%
  mutate(frag_id = paste0("TCs_", chr, "_", start, "_", end))

# =============================================================================
# 4. Annotation / strand count bar charts
# =============================================================================
 
annotation_counts <- data.frame(
  panel      = c("Strand", "Strand", "Annotation", "Annotation", "Annotation"),
  bar_color  = c("#fe0000", "#0725ff", "#565656", "#565656", "#565656"),
  annotation = c("Forward", "Reverse", "Promoters", "Enhancers", "ATACseq peaks"),
  count      = c(
    nrow(filter(fragment_annotated, strand == "+")),
    nrow(filter(fragment_annotated, strand == "-")),
    nrow(filter(fragment_annotated, islets_promoter == 1)),
    nrow(filter(fragment_annotated, islets_enhancer == 1)),
    nrow(filter(fragment_annotated, islet_atac == 1))
  )
)
annotation_counts$annotation <- factor(
  annotation_counts$annotation,
  levels = c("Forward", "Reverse", "Promoters", "Enhancers", "ATACseq peaks")
)
 
# Shared theme for both panels
bar_theme <- theme_minimal(base_size = 14) +
  theme(
    axis.text.x  = element_text(size = 10, color = "black"),
    axis.title.x = element_text(size = 14, color = "black")
  )
 
strand_bar_plot <- ggplot(
  filter(annotation_counts, panel == "Strand"),
  aes(x = annotation, y = count, fill = bar_color)
) +
  geom_col(color = "black") +
  scale_fill_identity() +
  bar_theme +
  theme(
    axis.text.y  = element_text(size = 12, color = "black"),
    axis.title.y = element_text(size = 12, color = "black")
  ) +
  labs(x = "CAGE Strandedness", y = "Number of Fragments") +
  ylim(0, 1200)
 
chromhmm_bar_plot <- ggplot(
  filter(annotation_counts, panel == "Annotation"),
  aes(x = annotation, y = count, fill = bar_color)
) +
  geom_col(color = "black") +
  scale_fill_identity() +
  bar_theme +
  theme(axis.text.y = element_blank()) +
  labs(x = "ChromHMM Annotations", y = "") +
  ylim(0, 1200)
 
counts_combined_plot <- plot_grid(
  strand_bar_plot, NULL, chromhmm_bar_plot,
  align       = "h",
  axis        = "l",
  rel_widths  = c(1.5, 0, 2),
  ncol        = 3,
  nrow        = 1
)
ggsave(
  plot     = counts_combined_plot,
  file     = "figures/annotation_strand_count.png",
  dpi      = 300,
  width    = 6.5,
  height   = 3
)
 
 # =============================================================================
# 5. Raincloud plot: z-score by strand concordance × element type
# =============================================================================
 
# Label fragments by strand concordance (concordant = "+" strand, discordant = "-")
raincloud_df <- fragment_annotated %>%
  mutate(
    strand_concordance = ifelse(
      strand == "+", "CAGE-strand concordant", "CAGE-strand discordant"
    ),
    element_type = case_when(
      islets_promoter == 1 ~ "promoters",
      islets_enhancer == 1 ~ "enhancers",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(element_type)) %>%
  mutate(
    category = paste(strand_concordance, element_type),
    category = factor(category, levels = c(
      "CAGE-strand concordant promoters",
      "CAGE-strand concordant enhancers",
      "CAGE-strand discordant promoters",
      "CAGE-strand discordant enhancers"
    ))
  )
 
# Colors map concordance direction (red = concordant, blue = discordant)
category_colors <- c(
  "CAGE-strand concordant promoters" = "red",
  "CAGE-strand concordant enhancers" = "red",
  "CAGE-strand discordant promoters" = "blue",
  "CAGE-strand discordant enhancers" = "blue"
)
 
raincloud_plot <- ggplot(
  raincloud_df,
  aes(x = category, y = zscore, fill = category, color = category)
) +
  # Half-violin density (right side); slab spans 25–75% quantile interval
  stat_halfeye(
    adjust       = 0.7,
    width        = 0.75,
    .width       = 0.5,
    point_colour = NA,
    interval_colour = NA,
    slab_alpha   = 0.6,
    position     = position_nudge(x = 0.1)
  ) +
  # Jittered raw data points (left side)
  geom_jitter(
    size     = 1,
    alpha    = 0.02,
    position = position_jitternudge(
      width       = 0.1,
      height      = 0,
      x           = -0.2,
      nudge.from  = "jittered.x"
    )
  ) +
  # Thin box plot for reference (centre)
  geom_boxplot(
    width        = 0.08,
    outlier.shape = NA,
    alpha        = 0.5,
    color        = "grey30",
    fill         = "white"
  ) +
  scale_fill_manual(values  = category_colors) +
  scale_color_manual(values = category_colors) +
  labs(x = NULL, y = "Fragment activity Z-score") +
  theme_minimal(base_size = 12) +
  theme(
    legend.position   = "none",
    axis.text.x       = element_blank(),
    panel.background  = element_rect(fill = "transparent", color = NA),
    plot.background   = element_rect(fill = "transparent", color = NA),
    plot.title        = element_text(face = "bold", size = 14, margin = margin(b = 10))
  )
 
ggsave(
  plot     = raincloud_plot,
  filename = "figures/density_plot.png",
  dpi      = 300,
  width    = 3.44,
  height   = 2.56
)

# =============================================================================
# 6. Univariate GLMs predicting fragment z-score from individual annotations
# =============================================================================
 
#' Fit a univariate Gaussian GLM and return a summary data frame.
#'
#' @param predictor  Character string naming the predictor column in `data`.
#' @param data       Data frame containing `zscore` and `predictor`.
#' @return           Data frame with columns: Factor, Coef, Std.Error, p-value,
#'                   t-value, Predictor, R2.
run_univariate_glm <- function(predictor, data) {
  fit <- glm(zscore ~ get(predictor, data), data = data)
 
  glm_summary <- coef(summary(fit))
 
  result <- data.frame(
    Factor    = rownames(glm_summary),
    Coef      = glm_summary[, 1],
    Std.Error = glm_summary[, 2],
    `t-value` = glm_summary[, 3],
    `p-value` = glm_summary[, 4],
    check.names = FALSE
  ) %>%
    mutate(
      Predictor = predictor,
      R2        = 1 - fit$deviance / fit$null.deviance
    )
 
  return(result)
}
 
# Predictor variables (annotation columns) to test
annotation_predictors <- c("strand", "islets_promoter", "islets_enhancer", "islet_atac")
 
# Fit GLMs and combine results; drop the intercept rows
glm_results <- lapply(annotation_predictors, run_univariate_glm, fragment_annotated) %>%
  do.call(rbind, .) %>%
  filter(Factor != "(Intercept)")
 
# Multiple-testing correction and significance labels
glm_results <- glm_results %>%
  mutate(
    # Benjamini-Hochberg FDR correction (primary)
    p.adj.BH         = p.adjust(`p-value`, method = "BH"),
    # Bonferroni correction (more conservative; shown for reference)
    p.adj.bonferroni = p.adjust(`p-value`, method = "bonferroni"),
    Significance     = case_when(
      p.adj.BH < 0.001 ~ "***",
      p.adj.BH < 0.01  ~ "**",
      p.adj.BH < 0.05  ~ "*",
      p.adj.BH < 0.1   ~ ".",
      TRUE             ~ ""
    )
  )
 
# Human-readable predictor labels and display order
predictor_labels <- c(
  "islets_promoter" = "Islet promoter",
  "islets_enhancer" = "Islet enhancer",
  "islet_atac"      = "Islet ATAC",
  "strand"          = "Forward CAGE strand"
)
glm_results <- glm_results %>%
  mutate(
    Predictor = recode(Predictor, !!!predictor_labels),
    color     = case_when(
      Predictor == "Forward CAGE strand" ~ "#fe0000",
      Predictor == "Islet promoter"      ~ "#9d4cd8",
      Predictor == "Islet enhancer"      ~ "#ffa400",
      Predictor == "Islet ATAC"          ~ "#565656"
    ),
    Predictor = factor(Predictor, levels = c(
      "Islet ATAC", "Islet enhancer", "Islet promoter", "Forward CAGE strand"
    ))
  )
 
# Forest-style dot plot of GLM t-values
glm_tvalue_plot <- ggplot(glm_results, aes(x = `t-value`, y = Predictor)) +
  geom_vline(aes(xintercept = 0)) +
  geom_col(aes(fill = color), width = 0.1) +
  geom_point(aes(color = color), size = 8) +
  # Significance labels: right of bar for positive t-values, left for negative
  geom_text(
    data    = filter(glm_results, `t-value` > 0),
    mapping = aes(label = Significance),
    size    = 10, hjust = -0.4, vjust = 0.75
  ) +
  geom_text(
    data    = filter(glm_results, `t-value` < 0),
    mapping = aes(label = Significance),
    size    = 10, hjust = 1.6, vjust = 0.75, color = "black"
  ) +
  scale_fill_identity() +
  scale_color_identity() +
  xlim(-5, 10) +
  labs(x = "GLM t-value", y = "") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y      = element_text(color = "black", size = 14),
    axis.text.x      = element_text(color = "black", size = 14),
    axis.title.x     = element_text(color = "black", size = 16)
  )
 
ggsave(
  plot   = glm_tvalue_plot,
  file   = "figures/glm_results.png",
  dpi    = 300,
  width  = 6.5,
  height = 4
)

# =============================================================================
# 7. GLM model diagnostics
# =============================================================================
 
# Re-fit GLM objects (one per predictor) for diagnostic extraction
glm_objects <- setNames(
  lapply(annotation_predictors, function(p) {
    glm(zscore ~ get(p, fragment_annotated), data = fragment_annotated)
  }),
  annotation_predictors
)
 
# Identify high-residual strand-discordant fragments (> mean + 2 SD)
fragment_annotated <- fragment_annotated %>%
  mutate(
    strand_residual = residuals(glm_objects[["strand"]], type = "deviance"),
    strand_label    = factor(
      ifelse(strand == "+", "strand\nconcordant", "strand\ndiscordant"),
      levels = c("strand\nconcordant", "strand\ndiscordant")
    )
  )
 
residual_cutoff  <- mean(fragment_annotated$strand_residual, na.rm = TRUE) +
                    2 * sd(fragment_annotated$strand_residual, na.rm = TRUE)
 
high_residual_discordant <- fragment_annotated %>%
  filter(strand == "-", strand_residual > residual_cutoff)
 
# Box + jitter plot highlighting high-residual discordant fragments
outlier_plot <- ggplot(fragment_annotated, aes(x = strand_label, y = zscore)) +
  geom_boxplot() +
  geom_jitter(alpha = 0.2, width = 0.2) +
  geom_point(data = high_residual_discordant, color = "red") +
  labs(x = "") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(color = "black"),
    axis.text.y = element_text(color = "black")
  )
ggsave("figures/outlier_fragments.png", outlier_plot, width = 4, height = 3, dpi = 300)
 
#' Build a 4-panel diagnostic plot for a fitted GLM.
#'
#' @param predictor  Name of the predictor (used as subtitle label).
#' @param fit        Fitted glm object.
#' @return           A cowplot grid of 4 ggplot panels.
make_diagnostic_plots <- function(predictor, fit) {
 
  residuals_vec <- residuals(fit, type = "response")
  fitted_vals   <- fitted(fit)
  plot_df       <- na.omit(data.frame(fitted = fitted_vals, residual = residuals_vec))
  diag_theme    <- theme_minimal(base_size = 12) +
                   theme(panel.background = element_blank(),
                         panel.grid       = element_blank())
 
  # Panel 1: Residuals vs Fitted — checks linearity and homoscedasticity
  group_means_resid <- plot_df %>%
    group_by(fitted) %>%
    summarise(mean_resid = mean(residual), .groups = "drop")
 
  p1 <- ggplot(plot_df, aes(x = fitted, y = residual)) +
    geom_jitter(alpha = 0.2, size = 0.8, height = 0, width = 0.25) +
    geom_hline(yintercept = 0, colour = "red", linetype = "dashed") +
    geom_line(data = group_means_resid, aes(y = mean_resid),
              colour = "blue", linewidth = 0.8) +
    labs(title = "Residuals vs Fitted", subtitle = predictor,
         x = "Fitted values", y = "Residuals") +
    diag_theme
 
  # Panel 2: Normal Q-Q — checks normality of residuals
  p2 <- ggplot(plot_df, aes(sample = residual)) +
    stat_qq(alpha = 0.3, size = 0.8) +
    stat_qq_line(colour = "red") +
    labs(title = "Normal Q-Q", subtitle = predictor,
         x = "Theoretical quantiles", y = "Sample quantiles") +
    diag_theme
 
  # Panel 3: Scale-Location — checks homoscedasticity
  plot_df <- plot_df %>%
    mutate(sqrt_abs_std_resid = sqrt(abs(scale(residual)[, 1])))
 
  group_means_sl <- plot_df %>%
    group_by(fitted) %>%
    summarise(mean_sl = mean(sqrt_abs_std_resid), .groups = "drop")
 
  p3 <- ggplot(plot_df, aes(x = fitted, y = sqrt_abs_std_resid)) +
    geom_jitter(alpha = 0.2, size = 0.8, height = 0, width = 0.25) +
    geom_line(data = group_means_sl, aes(y = mean_sl),
              colour = "blue", linewidth = 0.8) +
    labs(title = "Scale-Location", subtitle = predictor,
         x = "Fitted values", y = expression(sqrt("|Std. residuals|"))) +
    diag_theme
 
  # Panel 4: Residual histogram with overlaid normal curve — visual normality check
  bw <- 2 * IQR(residuals_vec) / length(residuals_vec)^(1/3)  # Freedman-Diaconis
 
  p4 <- ggplot(plot_df, aes(x = residual)) +
    geom_histogram(aes(y = after_stat(density)),
                   binwidth = bw, fill = "steelblue", colour = "white", alpha = 0.7) +
    stat_function(fun  = dnorm,
                  args = list(mean = mean(residuals_vec), sd = sd(residuals_vec)),
                  colour = "red", linewidth = 0.9) +
    labs(title = "Residual distribution", subtitle = predictor,
         x = "Residuals", y = "Density") +
    diag_theme
 
  plot_grid(p1, p2, p3, p4, nrow = 2)
}
 
# Save one PNG per predictor
for (predictor in annotation_predictors) {
  diag_plot <- make_diagnostic_plots(predictor, glm_objects[[predictor]])
  out_path  <- glue("figures/glm_diagnostics_{predictor}.png")
  ggsave(out_path, diag_plot, width = 10, height = 8)
  cat(sprintf("Saved: %s\n", out_path))
}
 
# =============================================================================
# 8. Per-configuration GLMs
# =============================================================================
# Fragments were tested in 4 plasmid configurations (varying promoter and
# insert position). This section repeats the annotation GLMs on each
# configuration's activity scores.
 
config_files <- list.files("data/config_activities", full.names = TRUE)
 
#' Run annotation GLMs on a single plasmid configuration activity file and
#' return a forest-style t-value plot.
#'
#' @param filename  Path to the config activity TSV (columns: Frag, zscore).
#' @return          A ggplot object.
run_config_glms_and_plot <- function(filename) {
 
  # Parse promoter and position identifiers from filename
  name_parts <- unlist(str_split(filename, "_"))
  promoter   <- name_parts[9]
  position   <- unlist(str_split(name_parts[11], "\\."))[1]
 
  # Load config-specific z-scores and join annotation columns
  config_df <- suppressMessages(
    read_tsv(filename, skip = 1, col_names = c("frag_id", "zscore"))
  ) %>%
    left_join(select(fragment_annotated, -zscore), by = "frag_id")
 
  # Fit annotation GLMs for this configuration
  config_glm_results <- lapply(annotation_predictors, run_univariate_glm, config_df) %>%
    do.call(rbind, .) %>%
    filter(Factor != "(Intercept)")
 
  # Multiple-testing correction and significance labels
  config_glm_results <- config_glm_results %>%
    mutate(
      p.adj.BH     = p.adjust(`p-value`, method = "BH"),
      Significance = case_when(
        p.adj.BH < 0.001 ~ "***",
        p.adj.BH < 0.01  ~ "**",
        p.adj.BH < 0.05  ~ "*",
        p.adj.BH < 0.1   ~ ".",
        TRUE             ~ ""
      ),
      Predictor = recode(Predictor, !!!predictor_labels),
      color     = case_when(
        Predictor == "Forward CAGE strand" ~ "#fe0000",
        Predictor == "Islet promoter"      ~ "#9d4cd8",
        Predictor == "Islet enhancer"      ~ "#ffa400",
        Predictor == "Islet ATAC"          ~ "#565656"
      ),
      Predictor = factor(Predictor, levels = c(
        "Islet ATAC", "Islet enhancer", "Islet promoter", "Forward CAGE strand"
      ))
    )
 
  # Forest-style dot plot
  ggplot(config_glm_results, aes(x = `t-value`, y = Predictor)) +
    geom_vline(aes(xintercept = 0)) +
    geom_col(aes(fill = color), width = 0.1) +
    geom_point(aes(color = color), size = 5) +
    geom_text(
      data    = filter(config_glm_results, `t-value` > 0),
      mapping = aes(label = Significance),
      size    = 3, hjust = -0.4, vjust = 0.75
    ) +
    geom_text(
      data    = filter(config_glm_results, `t-value` < 0),
      mapping = aes(label = Significance),
      size    = 3, hjust = 2, vjust = 0.75, color = "black"
    ) +
    scale_fill_identity() +
    scale_color_identity() +
    xlim(-8, 10) +
    labs(x = "GLM t-value", y = "", title = glue("{promoter}-{position}")) +
    theme_minimal() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.y      = element_text(color = "black", size = 8),
      axis.text.x      = element_text(color = "black", size = 10),
      axis.title.x     = element_text(color = "black", size = 12)
    )
}
 
config_glm_plots <- lapply(config_files, run_config_glms_and_plot)
 
# Assemble panel title, subtitle, and per-config plots into one figure
config_title    <- ggdraw() +
  draw_label("GLM results predicting configuration-specific fragment activity",
             fontface = "bold", size = 10)
 
config_subtitle <- ggdraw() +
  draw_label("(** = p < 0.01, *** = p < 0.001)", size = 8, color = "black", hjust = 1.2)
 
config_combined_plot <- plot_grid(
  config_title,
  config_subtitle,
  plot_grid(plotlist = config_glm_plots, nrow = 2),
  ncol        = 1,
  rel_heights = c(0.08, 0.02, 1)
)
 
ggsave(
  plot     = config_combined_plot,
  filename = "figures/glm_results_per_config.png",
  dpi      = 300,
  width    = 6,
  height   = 6
)