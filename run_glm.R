# Load libraries
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(glmnet))
suppressPackageStartupMessages(library(cowplot))
set.seed(2023)

# Read in results
df <- read.table("./cage_quant/quant.tsv", header = TRUE)
df$coord <- rownames(df)
rownames(df) <- NULL
df <- separate_wider_delim(df, cols = coord, delim="_", names=c("TC", "chr", "start", "end"))
zscores.with.coords <- select(df , c(zscore, chr, start, end))
zscores.with.coords$start <- as.numeric(zscores.with.coords$start)
zscores.with.coords$end <- as.numeric(zscores.with.coords$end)

# Read in islet tissue chromHMM annotations
tissue.overlap.filenames <- list.files(path = "./data/2022_05_20_cage_tissue_overlap", pattern = "*_loj198bp.tsv", full.names = TRUE)
tissue.overlap.filenames <- tissue.overlap.filenames[grepl("slet", tissue.overlap.filenames)]

# Read in strand of origin annotations
strand.df <- suppressMessages(read_tsv("./data/oligo_coords_with_strand.tsv", col_names = c("chr", "start", "end", "chr2", "start2", "end2", "sequence", "strand"))) %>%
    select(chr, start, end, strand, sequence)
strand.df <- left_join(zscores.with.coords, strand.df)

zscores.with.coords.copy <- zscores.with.coords

# Create annotations for each tissue file for each insert. 
new_cols <- lapply(tissue.overlap.filenames, function(file){
    
    # Get annotation type from filename
    file.parts <- str_split(basename(file), "_")[[1]]
    tissue <- file.parts[3]
    annot <- file.parts[4]
    
    # Select location of oligo and number of complete overlaps with tissue annotations
    loj <- suppressMessages(read_tsv(file, col_names = c("chr", "start", "end", "count"))) %>% 
        select(c(chr, start, end, count))

    # Join tissue annotation file to zscore df by oligo location
    zscores.with.coords.copy <- inner_join(zscores.with.coords.copy, loj, by = c("chr" = "chr", "start" = "start", "end" = "end"))
    
    # Recode count column. If an overlap with a tissue annotation exists, value
    # is 1. Otherwise, value is 0, indicating no overlap.
    zscores.with.coords.copy$chr2 <- ifelse(zscores.with.coords.copy$count > 0, 1, 0)
    
    
    # Create new column name containing annotation and tissue type
    if (annot == "promoters"){
        new_col <- tolower(paste(tissue, "Promoter", sep = "_"))
    }
    else if (annot == "enhancers"){
        new_col <- tolower(paste(tissue, "Enhancer", sep = "_"))
    }
    else {
        new_col <- tolower(paste(tissue, "atac", sep = "_"))
    }
    # Name the annotation column
    names(zscores.with.coords.copy)[names(zscores.with.coords.copy) == "chr2"] <- new_col
    return(select(zscores.with.coords.copy, c(new_col))) 
})

# Bind all of the new annotation columns together
tissue_annots_df <- do.call("cbind", new_cols)

tissue.strands.annots <- cbind(strand.df, tissue_annots_df)
all.df <- select(tissue.strands.annots, -c(sequence))
all.df <- left_join(all.df, zscores.with.coords)

# Run GLM models
run.single.models <- function(var){
    df_reg <- glm(zscore ~ get(var, all.df), data = all.df)
    df_reg_std <- rownames_to_column(as.data.frame((coef(summary(df_reg))[,2])))
    df_reg_t <- rownames_to_column(as.data.frame((coef(summary(df_reg))[,3])))
    df_reg_p <- rownames_to_column(as.data.frame((coef(summary(df_reg))[,4])))
    df_reg_coefs <- rownames_to_column(as.data.frame((coef(df_reg))))
    colnames(df_reg_coefs) <- c("Factor", "Coef")
    colnames(df_reg_std) <- c("Factor", "Std.Error")
    colnames(df_reg_t) <- c("Factor", "t-value")
    colnames(df_reg_p) <- c("Factor", "p-value")
    df_reg_coefs <- suppressMessages(left_join(df_reg_coefs, df_reg_std))
    df_reg_coefs <- suppressMessages(left_join(df_reg_coefs, df_reg_p))
    df_reg_coefs <- suppressMessages(left_join(df_reg_coefs, df_reg_t))
    df_reg_coefs <- mutate(df_reg_coefs, Predictor = var)
    return(df_reg_coefs)
}

vars <- c("strand", "islets_promoter", "islets_enhancer", "islet_atac")

models <- lapply(vars, run.single.models)
models <- do.call("rbind", models) %>% subset(Factor != "(Intercept)") %>% select(-Factor) %>% mutate("Significance" = ifelse(`p-value` < 0.001, "***", 
                                                                                                                      ifelse(`p-value` < 0.01, "**",
                                                                                                                      ifelse(`p-value` < 0.05, "*",
                                                                                                                      ifelse(`p-value` < 0.1, ".", "")))))
