# Load Libraries
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(glmnet))
set.seed(2022)

# Load motifs
subset.motifs <- scan("./cage_wald/subset_best_motifs/motif_list.encode_1KG.trimmed.txt", character())
subset.motifs <- c('zscore', subset.motifs)

all.motifs <- scan("./cage_wald/subset_best_motifs/trimmed_clusters_encode_motifs_1000G_list.csv", character())
all.motifs <- c('zscore', all.motifs)

chrom_states <- c('8_Genic_enhancer', 'X9_Active_enhancer_1', 'stretchEnhancer', 'typicalEnhancer',
                  'X10_Active_enhancer_2', 'X11_Weak_enhancer', 'X14_Bivalent_poised_TSS', 'X16_Repressed_polycomb',
                  'X17_Weak_repressed_polycomb', 'X18_Quiescent_low_signal', 'X1_Active_TSS', 'X2_Weak_TSS',
                  'X3_Flanking_TSS', 'X5_Strong_transcription', 'X6_Weak_transcription', 'X8_Genic_enhancer',
                  'X9_Active_enhancer', 'X9_Active_enhancer_1', 'stretchEnhancer', 'typicalEnhancer', 'atac_seq')

# Promoter bias model
prom.df <- read_tsv("./cage_wald/wald_prom.zscores.tsv")
prom.df.subset <- select(prom.df, which(names(prom.df) %in% subset.motifs))

# Select x and y for promoter model
prom.x <- model.matrix(zscore ~ . , prom.df.subset)[,-1]
prom.y <- prom.df.subset$zscore

n <- nrow(prom.df)
train_rows <- sample(1:n, 0.8 * n, replace = F)

prom.x.test <- prom.x[-train_rows, ]
prom.x.train <- prom.x[train_rows, ]

prom.y.test <- prom.y[-train_rows]
prom.y.train <- prom.y[train_rows]


prom.models <- list()
for (i in 0:20) {
    name <- paste0("alpha", i/20)
    prom.models[[name]] <-
        cv.glmnet(prom.x, prom.y, type.measure="mse", alpha = i/20, family="gaussian")
}

prom.results <- data.frame()
for (i in 0:20) {
  name <- paste0("alpha", i/20)
  
  ## Use each model to predict 'y' given the Testing dataset
  predicted <- predict(prom.models[[name]], 
                       s=prom.models[[name]]$lambda.min, newx=prom.x.test)
  
  ## Calculate the Mean Squared Error...
  mse <- mean((prom.y.test - predicted)^2)
  
  ## Store the results
  temp <- data.frame(alpha=i/20, mse=mse, name=name)
  prom.results <- rbind(prom.results, temp)
}

final_prom_model <- prom.results[which.min(prom.results$mse), ]$alpha
final_prom_model <- prom.models[[paste0("alpha", final_prom_model)]]

prom.final.results <- as.data.frame(as.matrix(coef(final_prom_model, s = "lambda.min")))
prom.final.results$motif <- rownames(prom.final.results)
rownames(prom.final.results) <- NULL
prom.final.results <- subset(prom.final.results, s1 != 0 & motif != "(Intercept)")

write_tsv(prom.final.results, "./prom.results.final.tsv")

# Position bias model
pos.df <- read_tsv("./cage_wald/wald_pos.zscores.tsv")
pos.df.subset <- select(pos.df, which(names(pos.df) %in% subset.motifs))

pos.x <- model.matrix(zscore ~ . , pos.df.subset)[,-1]
pos.y <- pos.df.subset$zscore

n <- nrow(pos.df)
train_rows <- sample(1:n, 0.7 * n, replace = F)

pos.x.test <- pos.x[-train_rows, ]
pos.x.train <- pos.x[train_rows, ]

pos.y.test <- pos.y[-train_rows]
pos.y.train <- pos.y[train_rows]

pos.models <- list()
for (i in 0:20) {
    name <- paste0("alpha", i/20)
    
    pos.models[[name]] <-
        cv.glmnet(pos.x, pos.y, type.measure="mse", alpha = i/20, family="gaussian")
}

pos.results <- data.frame()
for (i in 0:20) {
  name <- paste0("alpha", i/20)
  
  ## Use each model to predict 'y' given the Testing dataset
  predicted <- predict(pos.models[[name]], 
                       s=pos.models[[name]]$lambda.1se, newx=pos.x.test)
  
  ## Calculate the Mean Squared Error...
  mse <- mean((pos.y.test - predicted)^2)
  
  ## Store the results
  temp <- data.frame(alpha=i/20, mse=mse, name=name)
  pos.results <- rbind(pos.results, temp)
}

final_pos_model <- pos.results[which.min(pos.results$mse), ]$alpha
final_pos_model <- pos.models[[paste0("alpha", final_pos_model)]]

pos.final.results <- as.data.frame(as.matrix(coef(final_pos_model, s = "lambda.min")))
pos.final.results$motif <- rownames(pos.final.results)
rownames(pos.final.results) <- NULL
pos.final.results <- subset(pos.final.results, s1 != 0 & motif != "(Intercept)")

write_tsv(pos.final.results, "./pos.results.final.tsv")