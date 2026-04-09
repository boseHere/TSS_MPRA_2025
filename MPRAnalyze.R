#' Perform comparative MPRA analysis on a selected subset of features
library(tidyverse)
library(MPRAnalyze)
library(ggdist)
library(cowplot)

# Paths
INPUT_OBJECT <- "data/quant.rds"
OUTPUT_DIR   <- "results"

# Subsets to run
SUBSETS <- c("up", "down", "INS", "SCP1")

# Create output directory if needed
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# Load the MPRA object (already depth-corrected)
# =============================================================================
t        <- load(INPUT_OBJECT, verbose = TRUE)
full.obj  <- get(t)
obj       <- full.obj   

# =============================================================================
# Helper: run comparative analysis for one subset
# =============================================================================
run_subset <- function(obj, subset_name, out_dir) {

  cat("\n", strrep("=", 60), "\n")
  cat("  Subset:", subset_name, "\n")
  cat(strrep("=", 60), "\n")

  #' Filters the MPRA object to rows/columns matching `subset_name`,
  #' constructs a subset-specific `MpraObject`, and runs a comparative
  #' differential analysis using `analyzeComparative()`.
  curr_dnaCounts <- as.data.frame(obj@dnaCounts) %>%
    dplyr::select(contains(subset_name)) %>% as.matrix()
  curr_rnaCounts <- as.data.frame(obj@rnaCounts) %>%
    dplyr::select(contains(subset_name)) %>% as.matrix()
  curr_dnaAnnot  <- obj@dnaAnnot[grepl(subset_name, rownames(obj@dnaAnnot)), ]
  curr_rnaAnnot  <- obj@rnaAnnot[grepl(subset_name, rownames(obj@rnaAnnot)), ]
  curr_controls  <- obj@controls[grepl(subset_name, rownames(obj@dnaAnnot))]

  curr.obj <- MpraObject(
    dnaCounts = curr_dnaCounts,
    rnaCounts = curr_rnaCounts,
    dnaAnnot  = curr_dnaAnnot,
    rnaAnnot  = curr_rnaAnnot,
    controls  = curr_controls
  )

  #' The RNA design formula is chosen based on the subset
  if (subset_name %in% c("up", "down")) {
    model <- "~ rep + prom"
    curr_rnaAnnot$prom <- factor(curr_rnaAnnot$prom, levels = c("down", "up"))
  } else {
    model <- "~ rep + pos"
    curr_rnaAnnot$pos  <- factor(curr_rnaAnnot$pos,  levels = c("INS", "SCP1"))
  }

  curr.obj <- analyzeComparative(
    obj           = curr.obj,
    dnaDesign     = ~ barcode,
    rnaDesign     = as.formula(model),
    reducedDesign = ~ rep,
    fit.se        = TRUE
  )

  res         <- testLrt(curr.obj)
  res$Frag    <- rownames(res)
  res$subset  <- subset_name
  res$padj    <- p.adjust(res$pval, method = "fdr")
  res$sig     <- res$padj < 0.05
  rownames(res) <- NULL

  out_file <- file.path(out_dir, paste0("lrt_", subset_name, ".tsv"))
  write_tsv(as.data.frame(res), out_file)
  cat("  Results written to:", out_file, "\n")
  cat("  Significant fragments (FDR < 0.05):", sum(res$sig, na.rm = TRUE),
      "/", nrow(res), "\n")
}

# =============================================================================
# Run all four subsets
# =============================================================================
results_list <- lapply(SUBSETS, function(s) run_subset(obj, s, OUTPUT_DIR))
names(results_list) <- SUBSETS

# Combined results table
all_results <- bind_rows(results_list)
write_tsv(all_results, file.path(OUTPUT_DIR, "lrt_all_subsets.tsv"))
cat("\nCombined results written to:", file.path(OUTPUT_DIR, "lrt_all_subsets.tsv"), "\n")

# =============================================================================
# Global interaction test (full 2x2 model)
# =============================================================================
# Full RNA model: ~ promoter + position + promoter:position
# Reduced model (for interaction test): ~ promoter + position
# This directly tests whether the promoter effect DEPENDS on position.

cat("\n==== INTERACTION TEST: promoter x position ====\n")
obj_full <- analyzeComparative(
  obj           = obj,           # full dataset
  dnaDesign     = ~ barcode,
  rnaDesign     = ~ rep + prom * pos,   # full factorial
  reducedDesign = ~ rep + prom + pos,   # additive (no interaction)
  correctControls = TRUE,
  verbose = TRUE
)
res_interaction <- testLrt(obj_full)
res_interaction$fragment <- rownames(res_interaction)
cat("Fragments with significant INTERACTION (FDR < 0.05):",
    sum(res_interaction$fdr < 0.05, na.rm = TRUE), "/", nrow(res_interaction), "\n")
cat("logFC (interaction term) summary:\n")
print(summary(res_interaction$logFC))

# Also test promoter main effect (in full model context)
cat("\n==== PROMOTER MAIN EFFECT (additive model) ====\n")
obj_prom_main <- analyzeComparative(
  obj           = obj,
  dnaDesign     = ~ barcode,
  rnaDesign     = ~ rep + prom + pos,
  reducedDesign = ~ rep + pos,
  correctControls = TRUE,
  verbose = FALSE
)
res_prom_main <- testLrt(obj_prom_main)
res_prom_main$fragment <- rownames(res_prom_main)
cat("Significant promoter main effect (FDR < 0.05):",
    sum(res_prom_main$fdr < 0.05, na.rm = TRUE), "/", nrow(res_prom_main), "\n")

# Test position main effect (in full model context)
cat("\n==== POSITION MAIN EFFECT (additive model) ====\n")
obj_pos_main <- analyzeComparative(
  obj           = obj,
  dnaDesign     = ~ barcode,
  rnaDesign     = ~ rep + prom + pos,
  reducedDesign = ~ rep + prom,
  correctControls = TRUE,
  verbose = FALSE
)
res_pos_main <- testLrt(obj_pos_main)
res_pos_main$fragment <- rownames(res_pos_main)
cat("Significant position main effect (FDR < 0.05):",
    sum(res_pos_main$fdr < 0.05, na.rm = TRUE), "/", nrow(res_pos_main), "\n")

# Write results
write_tsv(as.data.frame(res_interaction), file.path(OUTPUT_DIR, "lrt_interaction.tsv"))
write_tsv(as.data.frame(res_prom_main),   file.path(OUTPUT_DIR, "lrt_prom_main.tsv"))
write_tsv(as.data.frame(res_pos_main),    file.path(OUTPUT_DIR, "lrt_pos_main.tsv"))
cat("\nResults written to:", OUTPUT_DIR, "\n")

# =============================================================================
# Subset model visualisation
# =============================================================================
# Use the combined per-subset results table built above
df <- all_results

# ---------------------------------------------------------------------------
# Plot 1 — Promoter bias within positional configurations (up / down)
# ---------------------------------------------------------------------------
pos.sub  <- subset(df, subset %in% c("up", "down"))
palette1 <- c("#bfbfbf", "#1a3b69")

p1 <- ggplot(pos.sub, aes(x = subset, y = logFC)) +

  # 1. Half-violin (density) — significant fragments only
  ggdist::stat_halfeye(
    data          = subset(pos.sub, sig),
    aes(x = subset, y = logFC),
    adjust        = 0.6,
    width         = 0.5,
    justification = -0.25,
    .width        = 0,
    point_colour  = NA,
    alpha         = 0.75,
    side          = "right",
    fill          = "#bdd3ef",
    color         = "#102541"
  ) +

  # 2. Jittered raw points
  geom_jitter(aes(color = sig),
    width  = 0.08,
    height = 0,
    size   = 1.5,
    alpha  = 0.45,
    shape  = 16
  ) +

  scale_fill_manual(values   = palette1) +
  scale_colour_manual(values = palette1) +

  labs(
    title = "Promoter bias across position\nsubsets",
    x     = "Position",
    y     = "logFC",
    color = "Significant (5% FDR)"
  ) +

  theme_classic(base_size = 14) +
  theme(
    legend.position    = "bottom",
    plot.title         = element_text(face = "bold", size = 13),
    plot.subtitle      = element_text(colour = "grey40", margin = margin(b = 10)),
    axis.text.x        = element_text(face = "bold"),
    panel.grid.major.y = element_line(colour = "grey90"),
    plot.margin        = margin(0, 0, 0, 0),
    plot.background    = element_blank(),
    panel.background   = element_blank()
  ) +
  guides(
    color = guide_legend(
      title.position = "top",
      override.aes   = list(alpha = 1)
    )
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.5)

print(p1)

# ---------------------------------------------------------------------------
# Plot 2 — Position bias within promoter configurations (INS / SCP1)
# ---------------------------------------------------------------------------
prom.sub <- subset(df, subset %in% c("INS", "SCP1"))
palette2 <- c("#bfbfbf", "#e04009")

p2 <- ggplot(prom.sub, aes(x = subset, y = logFC)) +

  # 1. Half-violin (density) — significant fragments only
  ggdist::stat_halfeye(
    data          = subset(prom.sub, sig),
    aes(x = subset, y = logFC),
    adjust        = 0.6,
    width         = 0.5,
    justification = -0.25,
    .width        = 0,
    point_colour  = NA,
    alpha         = 0.75,
    side          = "right",
    fill          = "#faa689",
    color         = "#b03207"
  ) +

  # 2. Jittered raw points
  geom_jitter(aes(color = sig),
    width  = 0.08,
    height = 0,
    size   = 1.5,
    alpha  = 0.45,
    shape  = 16
  ) +

  scale_fill_manual(values   = palette2) +
  scale_colour_manual(values = palette2) +

  labs(
    title = "Position bias across promoter\nsubsets",
    x     = "Promoter",
    y     = "logFC",
    color = "Significant (5% FDR)"
  ) +

  theme_classic(base_size = 14) +
  theme(
    legend.position    = "bottom",
    plot.title         = element_text(face = "bold", size = 13),
    plot.subtitle      = element_text(colour = "grey40", margin = margin(b = 10)),
    axis.text.x        = element_text(face = "bold"),
    panel.grid.major.y = element_line(colour = "grey90"),
    plot.margin        = margin(0, 0, 0, 0),
    plot.background    = element_blank(),
    panel.background   = element_blank()
  ) +
  guides(
    color = guide_legend(
      title.position = "top",
      override.aes   = list(alpha = 1)
    )
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.5)

print(p2)

# ---------------------------------------------------------------------------
# Combined panel
# ---------------------------------------------------------------------------
options(repr.plot.width = 8.267, repr.plot.height = 4.606)
final_plot <- plot_grid(
  plotlist   = list(NULL, p1, NULL, p2),
  nrow       = 1,
  rel_widths = c(0.1, 0.4, 0.1, 0.4)
)

final_plot
ggsave(
  plot     = final_plot,
  filename = file.path(OUTPUT_DIR, "subset_models.png"),
  dpi      = 300,
  width    = 8.267,
  height   = 4.606
)
cat("\nFigure written to:", file.path(OUTPUT_DIR, "subset_models.png"), "\n")

#' Output:
#' - Per-subset TSV files:       results/lrt_up.tsv, lrt_down.tsv, lrt_INS.tsv, lrt_SCP1.tsv
#' - Combined subsets TSV:       results/lrt_all_subsets.tsv
#' - Interaction test TSV:       results/lrt_interaction.tsv
#' - Promoter main effect TSV:   results/lrt_prom_main.tsv
#' - Position main effect TSV:   results/lrt_pos_main.tsv
#' - Combined figure (PNG):      results/subset_models.png