## Load libraries
library(dplyr)
library(stringr)
library(gProfileR)
library(Hmisc)
library(corrplot)
library(pheatmap)
library(dendsort)


## Set working directory
setwd("C:/Users/Tony Lin/Desktop/Wiita_lab/Computational_tools/random_r_scripts/arun-proteostais/")


## Read DepMap (19Q1) Expression File (57820 genes/1165 cell lines)
exp = read.csv("CCLE_depMap_19Q1_TPM.csv", stringsAsFactors = F)


## Read DepMap (19Q1) CRISPR File (17634 genes/558 cell lines)
crispr = read.csv("gene_effect_corrected_19Q1.csv", stringsAsFactors = F)


## Read proteasome genes
proteasome = read.csv("all-proteasome-genes.csv", stringsAsFactors = F)[[1]]
proteasome = toupper(proteasome)
proteasome = unique(proteasome)
write.csv(proteasome, "unique-proteostasis-genes.csv")


## Identify human MM cell lines
meta = read.csv("DepMap-2019q1-celllines_v2.csv", stringsAsFactors = F) 
MM_meta = meta %>%
  filter(Primary.Disease == "Myeloma") %>%
  mutate(HMCL = sub("_.*", "", CCLE_Name))
#save.image("proteostasis.RData")


## Clean up CRISPR screen data
#load("proteostasis.RData")
# Transpose crispr data frame to obtain GENES vs CELL LINES
rownames(crispr) = crispr$X
crispr = select(crispr, -1)
crispr = t(crispr)
crispr = as.data.frame(crispr)

# Clean up gene labels
crispr$GENE = sub("\\.\\..*", "", rownames(crispr))


## Clean up expression data
# Transpose exp data frame to obtain GENES vs CELL LINES
rownames(exp) = exp$X
exp = select(exp, -1)
exp = t(exp)
exp = as.data.frame(exp)

# Clean up gene labels
exp$GENE = sub("\\.\\..*", "", rownames(exp))
exp$ENSG = str_extract(rownames(exp), "\\.\\..*")
exp$ENSG = sub("^\\.\\.", "", exp$ENSG)
exp$ENSG = sub("\\.$", "", exp$ENSG)


## Filter for genes involved in proteostasis in expression data
# All genes in our list are found in expression data
sum(proteasome %in% exp$GENE)

# PSMA1 and PSMA2 are duplicated (paralogs) in expression data
table(exp$GENE[exp$GENE %in% proteasome])
exp$ENSG[exp$GENE == "PSMA1"]   # ENSG00000256206 is a paralog of PSMA1
exp$ENSG[exp$GENE == "PSMA2"]   # ENSG00000256646 is a paralog of PSMA2

hist(as.numeric(exp[which(exp$ENSG == "ENSG00000129084"), 1:1165]))   # PSMA1 expression
hist(as.numeric(exp[which(exp$ENSG == "ENSG00000256206"), 1:1165]))   # PSMA1-paralog

hist(as.numeric(exp[which(exp$ENSG == "ENSG00000106588"), 1:1165]))   # PSMA2 expression
hist(as.numeric(exp[which(exp$ENSG == "ENSG00000256646"), 1:1165]))   # PSMA2-paralog expression

# Filter for proteostasis genes and PSMA1 (ENSG00000129084) and PSMA2 (ENSG00000106588)
proExp = exp %>%
  filter(GENE %in% proteasome) %>%
  filter(!(ENSG %in% c("ENSG00000256206", "ENSG00000256646")))
rownames(proExp) = proExp$GENE
proExp = proExp[1:1165]


## Filter for genes involved in proteostasis in CRISPR data
# 406/441 genes are found in CRISPR screen data
sum(proteasome %in% crispr$GENE)

# All 406 genes are present only once in CRISPR screen data
sum(crispr$GENE %in% proteasome)

# Filter for proteostasis genes in CRISPR data
proCrispr = crispr %>%
  filter(GENE %in% proteasome)
rownames(proCrispr) = proCrispr$GENE
proCrispr = proCrispr[1:558]


## Generate pairwise correlation plot on single data type
plot_corr = function(df, cor_type = "pearson", mark_insig = "n") {
  # df = data frame with genes in rows and cell lines in columns
  # cor_type = character of length 1 indicating the type of correlation (pearson or spearman)
  # mark_insig = "n" means ignore significance level and "pch" means cross out insignificant correlations
  cor_mat = rcorr(t(df), type = cor_type)
  
  corrplot(cor_mat$r, type = "full", method = "shade",
           order = "hclust", hclust.method = "complete", 
           p.mat = cor_mat$P, sig.level = 0.05, insig = mark_insig, 
           tl.col = "black", tl.srt = 45, tl.cex = 2, cl.cex = 36,
           col = colorRampPalette(c("blue", "white", "red"))(300))
}
pdf("crispr_pearson_correlogram_558lines.pdf", width = 200, height = 200)
plot_corr(proCrispr, cor_type = "pearson")
dev.off()

pdf("crispr_spearman_correlogram_558lines.pdf", width = 200, height = 200)
plot_corr(proCrispr, cor_type = "spearman")
dev.off()

pdf("RNA_pearson_correlogram_1165lines.pdf", width = 200, height = 200)
plot_corr(proExp, cor_type = "pearson")
dev.off()

pdf("RNA_spearman_correlogram_1165lines.pdf", width = 200, height = 200)
plot_corr(proExp, cor_type = "spearman")
dev.off()

pdf("crispr_pearson_correlogram_558lines_pvalue_0.05.pdf", width = 200, height = 200)
plot_corr(proCrispr, cor_type = "pearson", mark_insig = "pch")
dev.off()

pdf("crispr_spearman_correlogram_558lines_pvalue_0.05.pdf", width = 200, height = 200)
plot_corr(proCrispr, cor_type = "spearman", mark_insig = "pch")
dev.off()

pdf("RNA_pearson_correlogram_1165lines_pvalue_0.05.pdf", width = 200, height = 200)
plot_corr(proExp, cor_type = "pearson", mark_insig = "pch")
dev.off()

pdf("RNA_spearman_correlogram_1165lines_pvalue_0.05.pdf", width = 200, height = 200)
plot_corr(proExp, cor_type = "spearman", mark_insig = "pch")
dev.off()


## Myeloma-only
pdf("crispr_pearson_correlogram_16MMlines.pdf", width = 200, height = 200)
plot_corr(proCrispr[names(proCrispr) %in% MM_meta$DepMap_ID], 
          cor_type = "pearson")
dev.off()

pdf("crispr_spearman_correlogram_16MMlines.pdf", width = 200, height = 200)
plot_corr(proCrispr[names(proCrispr) %in% MM_meta$DepMap_ID], 
          cor_type = "spearman")
dev.off()

pdf("RNA_pearson_correlogram_27MMlines.pdf", width = 200, height = 200)
plot_corr(proExp[names(proExp) %in% MM_meta$DepMap_ID], 
          cor_type = "pearson")
dev.off()

pdf("RNA_spearman_correlogram_27MMlines.pdf", width = 200, height = 200)
plot_corr(proExp[names(proExp) %in% MM_meta$DepMap_ID], 
          cor_type = "spearman")
dev.off()

pdf("crispr_pearson_correlogram_16MMlines_pvalue_0.05.pdf", width = 200, height = 200)
plot_corr(proCrispr[names(proCrispr) %in% MM_meta$DepMap_ID], 
          cor_type = "pearson", mark_insig = "pch")
dev.off()

pdf("crispr_spearman_correlogram_16MMlines_pvalue_0.05.pdf", width = 200, height = 200)
plot_corr(proCrispr[names(proCrispr) %in% MM_meta$DepMap_ID], 
          cor_type = "spearman", mark_insig = "pch")
dev.off()

pdf("RNA_pearson_correlogram_27MMlines_pvalue_0.05.pdf", width = 200, height = 200)
plot_corr(proExp[names(proExp) %in% MM_meta$DepMap_ID], 
          cor_type = "pearson", mark_insig = "pch")
dev.off()

pdf("RNA_spearman_correlogram_27MMlines_pvalue_0.05.pdf", width = 200, height = 200)
plot_corr(proExp[names(proExp) %in% MM_meta$DepMap_ID], 
          cor_type = "spearman", mark_insig = "pch")
dev.off()


# Compare correlation matrix between MM lines and all lines
## Generate pairwise correlation plot on single data type
compare_corr = function(df, df2) {
  # df = data frame with genes in rows and cell lines in columns (all cells)
  # df2 = data frame with genes in rows and cell lines in columsn (MM cells only)
  all_matP = rcorr(t(df), type = "pearson")$r
  MM_matP = rcorr(t(df2), type = "pearson")$r
  all_matS = rcorr(t(df), type = "spearman")$r
  MM_matS = rcorr(t(df2), type = "spearman")$r
  
  # Flatten correlation data
  ut = upper.tri(all_matP)
  all_matP = data.frame(row = rownames(all_matP)[row(all_matP)[ut]],
                       col = colnames(all_matP)[col(all_matP)[ut]],
                       all_cor_pearson = all_matP[ut])
  all_matS = data.frame(row = rownames(all_matS)[row(all_matS)[ut]],
                        col = colnames(all_matS)[col(all_matS)[ut]],
                        all_cor_spearman = all_matS[ut])
  MM_matP = data.frame(row = rownames(MM_matP)[row(MM_matP)[ut]],
                      col = colnames(MM_matP)[col(MM_matP)[ut]],
                      MM_cor_pearson = MM_matP[ut])
  MM_matS = data.frame(row = rownames(MM_matS)[row(MM_matS)[ut]],
                      col = colnames(MM_matS)[col(MM_matS)[ut]],
                      MM_cor_spearman = MM_matS[ut])
  
  # Clean data
  all_matP$gene_pairs = paste0(all_matP$row, "-", all_matP$col)
  all_matP = select(all_matP, -(1:2))
  all_matS$gene_pairs = paste0(all_matS$row, "-", all_matS$col)
  all_matS = select(all_matS, -(1:2))
  MM_matP$gene_pairs = paste0(MM_matP$row, "-", MM_matP$col)
  MM_matP = select(MM_matP, -(1:2))
  MM_matS$gene_pairs = paste0(MM_matS$row, "-", MM_matS$col)
  MM_matS = select(MM_matS, -(1:2))
  
  # Merge data frame
  merged = inner_join(all_matP, MM_matP, by = "gene_pairs") %>%
    inner_join(all_matS, by = "gene_pairs") %>%
    inner_join(MM_matS, by = "gene_pairs") %>%
    mutate(gene_pairs_rev = paste0(gene_pairs, "-", gene_pairs)) %>%
    mutate(gene_pairs_rev = str_extract(gene_pairs_rev, ".*-")) %>%
    mutate(gene_pairs_rev = str_extract(gene_pairs_rev, "-.*")) %>%
    mutate(gene_pairs_rev = sub("^-", "", gene_pairs_rev)) %>%
    mutate(gene_pairs_rev = sub("-$", "", gene_pairs_rev))
  
  return(merged)
}
corCrispr = compare_corr(proCrispr,
                         proCrispr[names(proCrispr) %in% MM_meta$DepMap_ID])
corExp = compare_corr(proExp,
                      proExp[names(proExp) %in% MM_meta$DepMap_ID])
write.csv(corCrispr, "crispr-correlation-values.csv", row.names = F)
write.csv(corExp, "exp-correlation-values.csv", row.names = F)


# Correlate genetic dependency vs gene expression
## Transpose data frames to obtain genes in columns and cell lines in rows
TproCrispr = as.data.frame(t(proCrispr))
TproExp = as.data.frame(t(proExp))

## Label columns as cirspr or exp data
colnames(TproCrispr) = paste0(colnames(TproCrispr), "_CRISPR")
colnames(TproExp) = paste0(colnames(TproExp), "_EXP")

## Merge data frames by cell lines
TproCrispr$cell_line = rownames(TproCrispr)
TproExp$cell_line = rownames(TproExp)
combo = inner_join(TproCrispr, TproExp, by = "cell_line")
rownames(combo) = combo$cell_line
combo = select(combo, -cell_line)
combo = as.data.frame(t(combo))

## Plot heatmaps
plot_corr2 = function(df, cor_type = "pearson", mark_insig = "n") {
  # df = data frame with genes in rows and cell lines in columns
  # cor_type = character of length 1 indicating the type of correlation (pearson or spearman)
  # mark_insig = "n" means ignore significance level and "pch" means cross out insignificant correlations
  cor_mat = rcorr(t(df), type = cor_type)
  cor_mat$r = cor_mat$r[grepl("_CRISPR", rownames(cor_mat$r)),
                        grepl("_EXP", colnames(cor_mat$r))]
  cor_mat$P = cor_mat$P[grepl("_CRISPR", rownames(cor_mat$P)),
                        grepl("_EXP", colnames(cor_mat$P))]

  # Filter out rows or columns with NA
  if (sum(apply(cor_mat$r, 1, function(x) all(is.nan(x)))) > 0) {
    filt = apply(apply(cor_mat$r, 1, function(x) all(is.nan(x))))
    cor_mat$r = cor_mat$r[!filt, ]
    cor_mat$P = cor_mat$P[!filt, ]
  }
  if (sum(apply(cor_mat$r, 2, function(x) all(is.nan(x)))) > 0) {
    filt = apply(cor_mat$r, 2, function(x) all(is.nan(x)))
    cor_mat$r = cor_mat$r[, !filt]
    cor_mat$P = cor_mat$P[, !filt]
  }
  
  rowSort = dist(cor_mat$r, method = "euclidean")   # Find row distance
  colSort = dist(t(cor_mat$r), method = "euclidean")   # Find column distance
  
  sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...)))  # For ordering dendrogram
  colors = colorRampPalette(c("darkblue", "white", "darkred"))(300)
  breakPTS = seq(from = -max(abs(cor_mat$r)),
                 to = max(abs(cor_mat$r)),
                 length.out = 300)
  
  # Unable to implement clustering using corrplot when matrix is asymmetrical
  pheatmap(cor_mat$r,
           cluster_cols = sort_hclust(hclust(colSort)),
           cluster_rows = sort_hclust(hclust(rowSort)),
           col = colors,
           display_numbers = round(cor_mat$P, 3),
           breaks = breakPTS,
           fontsize_row = 20,
           fontsize_col = 20
           )
}
pdf("crispr_vs_exp_pearson_heatmap_554lines.pdf", width = 200, height = 200)
plot_corr2(combo, cor_type = "pearson")
dev.off()

pdf("crisprvs_exp_spearman_heatmap_554lines.pdf", width = 200, height = 200)
plot_corr2(combo, cor_type = "spearman")
dev.off()

## Myeloma-only
pdf("crispr_vs_exp_pearson_heatmap_16MMlines.pdf", width = 200, height = 200)
plot_corr2(combo[names(combo) %in% MM_meta$DepMap_ID],
           cor_type = "pearson")
dev.off()

pdf("crispr_vs_exp_spearman_heatmap_16MMlines.pdf", width = 200, height = 200)
plot_corr2(combo[names(combo) %in% MM_meta$DepMap_ID],
           cor_type = "spearman")
dev.off()


# Compare correlation matrix between MM lines and all lines
## Generate pairwise correlation plot on single data type
compare_corr2 = function(df, df2) {
  # df = data frame with genes in rows and cell lines in columns (all cells)
  # df2 = data frame with genes in rows and cell lines in columsn (MM cells only)
  all_matP = rcorr(t(df), type = "pearson")$r
  all_matP = all_matP[grepl("_CRISPR", rownames(all_matP)),
                      grepl("_EXP", colnames(all_matP))]
  all_matS = rcorr(t(df), type = "spearman")$r
  all_matS = all_matS[grepl("_CRISPR", rownames(all_matS)),
                      grepl("_EXP", colnames(all_matS))]
  MM_matP = rcorr(t(df2), type = "pearson")$r
  MM_matP = MM_matP[grepl("_CRISPR", rownames(MM_matP)),
                    grepl("_EXP", colnames(MM_matP))]
  MM_matS = rcorr(t(df2), type = "spearman")$r
  MM_matS = MM_matS[grepl("_CRISPR", rownames(MM_matS)),
                    grepl("_EXP", colnames(MM_matS))]
  
  # Flatten correlation data
  ut = upper.tri(all_matP, diag = T)
  all_matP = data.frame(row = rownames(all_matP)[row(all_matP)[ut]],
                        col = colnames(all_matP)[col(all_matP)[ut]],
                        all_cor_pearson = all_matP[ut])
  all_matS = data.frame(row = rownames(all_matS)[row(all_matS)[ut]],
                        col = colnames(all_matS)[col(all_matS)[ut]],
                        all_cor_spearman = all_matS[ut])
  MM_matP = data.frame(row = rownames(MM_matP)[row(MM_matP)[ut]],
                       col = colnames(MM_matP)[col(MM_matP)[ut]],
                       MM_cor_pearson = MM_matP[ut])
  MM_matS = data.frame(row = rownames(MM_matS)[row(MM_matS)[ut]],
                       col = colnames(MM_matS)[col(MM_matS)[ut]],
                       MM_cor_spearman = MM_matS[ut])
  
  # Clean data
  all_matP$gene_pairs = paste0(all_matP$row, "-", all_matP$col)
  all_matP = select(all_matP, -(1:2))
  all_matS$gene_pairs = paste0(all_matS$row, "-", all_matS$col)
  all_matS = select(all_matS, -(1:2))
  MM_matP$gene_pairs = paste0(MM_matP$row, "-", MM_matP$col)
  MM_matP = select(MM_matP, -(1:2))
  MM_matS$gene_pairs = paste0(MM_matS$row, "-", MM_matS$col)
  MM_matS = select(MM_matS, -(1:2))
  
  # Merge data frame
  merged = inner_join(all_matP, MM_matP, by = "gene_pairs") %>%
    inner_join(all_matS, by = "gene_pairs") %>%
    inner_join(MM_matS, by = "gene_pairs")
  
  return(merged)
}
corCombo = compare_corr2(combo, combo[names(combo) %in% MM_meta$DepMap_ID])
write.csv(corCombo, "crispr-vs-exp-correlation-values.csv", row.names = F)


# Save data
#save.image("proteostasis-processed.RData")


