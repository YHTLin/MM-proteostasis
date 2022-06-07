### Survival Analysis
## Proteasome gene set (proteasome value)
load("C:/Users/Tony Lin/Desktop/Wiita_lab/Computational_tools/random_r_scripts/arun-proteostais/proteostasis-processed.RData")

## Set directory
setwd("C:/Users/Tony Lin/Desktop/Wiita_lab/Computational_tools/random_r_scripts/arun-commpass-analysis/")

## Read Salmon TPM
exp = read.delim("MMRF_CoMMpass_IA14a_E74GTF_Salmon_Gene_TPM.txt", stringsAsFactors = F)

### Survival Analysis for HSP proteins
## Convert ENSG to gene name
library(gProfileR)
convert_df = gconvert(exp$GENE_ID, organism = "hsapiens", target = "HGNC",
                      mthreshold = 1, filter_na = F)
exp$GENE = convert_df$target
exp$GENE_NAME = convert_df$description
exp$GENE_NAME = sub("\\[.*", "", exp$GENE_NAME)

## Clean up expression data
library(dplyr)
id_ND = names(exp)[grepl("MMRF_[0-9]{4}_1.*", names(exp))]   # Select newly diagnosed patients
exp_surv = exp[exp$GENE %in% proteasome, 
               c(id_ND, "GENE")]
row.names(exp_surv) = exp_surv$GENE
exp_surv = exp_surv[names(exp_surv) != "GENE"]
exp_surv = exp_surv[apply(exp_surv, 1, function(x) all(x > 1)), ]   # Remove genes with TPM < 1
exp_surv = as.data.frame(t(exp_surv))
exp_surv$public_id = sub("_._..$", "", row.names(exp_surv))
exp_surv = exp_surv %>%
  group_by(public_id) %>%
  summarise_each(mean)   # Average any duplicates (10 total)

## Read patient survival data
survival = read.csv("CoMMpass_IA14_FlatFiles/MMRF_CoMMpass_IA14_STAND_ALONE_SURVIVAL.csv", stringsAsFactors = F)
patientDF = read.csv("CoMMpass_IA14_FlatFiles/MMRF_CoMMpass_IA14_PER_PATIENT.csv", stringsAsFactors = F)

setwd("C:/Users/Tony Lin/Desktop/Wiita_lab/Computational_tools/random_r_scripts/arun-commpass-analysis/survival-salmon-proteostasis/survival-OS-salmon-proteostasis/")
for (i in names(exp_surv)[-1]) {
  ## Filter for patients with survival data
  survival2 = survival %>%
    dplyr::select(public_id, censos, oscdy) %>%   # Overall survival censor flag (censos) and date (oscdy)
    dplyr::rename(status = censos, time = oscdy) %>%   # Rename censos to status and oscdy to time
    #dplyr::select(public_id, censpfs, pfscdy) %>%   # Progression free survival censor flag (pfscdy) and date (oscdy)
    #dplyr::rename(status = censpfs, time = pfscdy) %>%   # Rename censos to status and oscdy to time      dplyr::filter(public_id %in% patientDF$PUBLIC_ID[patientDF$D_PT_PRIMARYREASON %in% c("Death", "")]) %>%   # Remove patients who exited program prematurely
    dplyr::inner_join(exp_surv, by = "public_id")
  
  ## Stratify by gene expression
  survival2$isTop = NA
  survival2$order = rank(survival2[i]) / nrow(survival2)
  survival2$isTop[survival2$order > 0.8] = T
  survival2$isTop[survival2$order < 0.2] = F
  
  ## Plot survival curve
  library(survival)
  library(survminer)
  fit <- survfit(Surv(time, status) ~ isTop, data = survival2)
  cairo_pdf(paste0(i, "-OS-salmon.pdf"), width = 6, height = 6)
  g <- ggsurvplot(fit, 
                  data = survival2, 
                  size = 0.5,                 # change line size
                  palette = 
                    c("#E7B800", "#2E9FDF"),# custom color palettes
                  conf.int = FALSE,          # Add confidence interval
                  pval = TRUE,              # Add p-value
                  pval.size = 6,
                  risk.table = T,        # Add risk table
                  risk.table.col = "strata",# Risk table color by groups
                  legend.labs = 
                    c("Bottom 20%", "Top 20%"),    # Change legend labels
                  #risk.table.height = 0.25, # Useful to change when you have multiple groups
                  ggtheme = theme_classic() + 
                    theme(text = element_text(size = 15),
                          axis.title = element_text(size = 18)),      # Change ggplot2 theme
                  xlab = "Days",
                  ylab = "Overall Survival",
                  title = i
  )   # default pdf with 6" x 6"
  print(g)
  dev.off()
}

output = data.frame(Gene = names(exp_surv)[-1],
                    median_OS_high_exp = rep(NA, ncol(exp_surv) - 1),
                    median_OS_low_exp = rep(NA, ncol(exp_surv) - 1),
                    median_OS_diff = rep(NA, ncol(exp_surv) - 1),
                    OS_pval = rep(NA, ncol(exp_surv) - 1),
                    median_PFS_high_exp = rep(NA, ncol(exp_surv) - 1),
                    median_PFS_low_exp = rep(NA, ncol(exp_surv) - 1),
                    median_PFS_diff = rep(NA, ncol(exp_surv) - 1),
                    PFS_pval = rep(NA, ncol(exp_surv) - 1))

# Flat file with all pvalues
for (i in names(exp_surv)[-1]) {
  ## Filter for patients with survival data
  survival2 = survival %>%
    dplyr::select(public_id, censos, oscdy) %>%   # Overall survival censor flag (censos) and date (oscdy)
    dplyr::rename(status = censos, time = oscdy) %>%   # Rename censos to status and oscdy to time
    #dplyr::select(public_id, censpfs, pfscdy) %>%   # Progression free survival censor flag (pfscdy) and date (oscdy)
    #dplyr::rename(status = censpfs, time = pfscdy) %>%   # Rename censos to status and oscdy to time      dplyr::filter(public_id %in% patientDF$PUBLIC_ID[patientDF$D_PT_PRIMARYREASON %in% c("Death", "")]) %>%   # Remove patients who exited program prematurely
    dplyr::inner_join(exp_surv, by = "public_id")
  
  ## Stratify by gene expression
  survival2$isTop = NA
  survival2$order = rank(survival2[i]) / nrow(survival2)
  survival2$isTop[survival2$order > 0.8] = T
  survival2$isTop[survival2$order < 0.2] = F
  
  ## Plot survival curve
  library(survival)
  library(survminer)
  fit <- survfit(Surv(time, status) ~ isTop, data = survival2)
  output$median_OS_high_exp[which(output$Gene == i)] <- surv_median(fit)[2, 2]
  output$median_OS_low_exp[which(output$Gene == i)] <- surv_median(fit)[1, 2]
  output$median_OS_diff[which(output$Gene == i)] <- surv_median(fit)[2, 2] - surv_median(fit)[1, 2]
  output$OS_pval[which(output$Gene == i)] <-surv_pvalue(fit)$pval
}

for (i in names(exp_surv)[-1]) {
  ## Filter for patients with survival data
  survival2 = survival %>%
    #dplyr::select(public_id, censos, oscdy) %>%   # Overall survival censor flag (censos) and date (oscdy)
    #dplyr::rename(status = censos, time = oscdy) %>%   # Rename censos to status and oscdy to time
    dplyr::select(public_id, censpfs, pfscdy) %>%   # Progression free survival censor flag (pfscdy) and date (oscdy)
    dplyr::rename(status = censpfs, time = pfscdy) %>%   # Rename censos to status and oscdy to time      dplyr::filter(public_id %in% patientDF$PUBLIC_ID[patientDF$D_PT_PRIMARYREASON %in% c("Death", "")]) %>%   # Remove patients who exited program prematurely
    dplyr::inner_join(exp_surv, by = "public_id")
  
  ## Stratify by gene expression
  survival2$isTop = NA
  survival2$order = rank(survival2[i]) / nrow(survival2)
  survival2$isTop[survival2$order > 0.8] = T
  survival2$isTop[survival2$order < 0.2] = F
  
  ## Plot survival curve
  library(survival)
  library(survminer)
  fit <- survfit(Surv(time, status) ~ isTop, data = survival2)
  output$median_PFS_high_exp[which(output$Gene == i)] <- surv_median(fit)[2, 2]
  output$median_PFS_low_exp[which(output$Gene == i)] <- surv_median(fit)[1, 2]
  output$median_PFS_diff[which(output$Gene == i)] <- surv_median(fit)[2, 2] - surv_median(fit)[1, 2]
  output$PFS_pval[which(output$Gene == i)] <-surv_pvalue(fit)$pval
}
write.csv(output, "survival-pvals-salmon-proteostasis.csv", row.names = F)



####################################################################
## Set directory
####################################################################
setwd("C:/Users/Tony Lin/Desktop/Wiita_lab/Computational_tools/random_r_scripts/arun-commpass-analysis/")

## Read pval file
x = read.csv("survival-pvals-salmon-proteostasis.csv", stringsAsFactors = F)
x$rank_OS = rank(x$OS_pval)
x$rank_PFS = rank(x$PFS_pval)

## Adjust for visualization
x$pval = x$OS_pval
x$rank = x$rank_OS

## Plot p-values
library(ggplot2)
library(ggrepel)
ggplot(x) +
  geom_point(aes(x = rank, y = -log10(pval)),
             alpha = 0.4) +
  geom_point(data = subset(x, grepl("PSMA", Gene)),
             aes(x = rank, y = -log10(pval)),
             color = "red", alpha = 0.6, size = 3) +
  geom_point(data = subset(x, grepl("PSMB", Gene)),
             aes(x = rank, y = -log10(pval)),
             color = "darkgreen", alpha = 0.6, size = 3) +
  geom_point(data = subset(x, Gene %in% c("HSPA9", "HSPA5", "HSPA4",
                                          "HSPA8", "HSPA13")),
             aes(x = rank, y = -log10(pval)),
             color = "blue", alpha = 0.6, size = 3) +
  geom_text_repel(
    data = subset(x, Gene %in% c("HSPA9", "HSPA5", "HSPA4",
                                 "HSPA8", "HSPA13")),
    aes(x = rank, y = -log10(pval), 
        label = Gene),
    color = "blue",
    size = 3, 
    direction = "x", 
    nudge_y = 0.5 + log10(subset(x, Gene %in% c("HSPA9", "HSPA5", "HSPA4",
                                                "HSPA8", "HSPA13"))$pval)) +
  geom_text_repel(
    data = subset(x, grepl("PSMA", Gene)),
    aes(x = rank, y = -log10(pval),
        label = Gene),
    color = "red",
    size = 3, 
    direction = "x", 
    nudge_y = 6.25 + log10(subset(x, grepl("PSMA", Gene))$pval)) +
  geom_text_repel(
    data = subset(x, grepl("PSMB", Gene)),
    aes(x = rank, y = -log10(pval),
        label = Gene),
    color = "darkgreen",
    size = 3, 
    direction = "y", 
    nudge_x = 230 - subset(x, grepl("PSMB", Gene))$rank) +
  xlab("P-value Rank") +
  ylab(expression("- log"[10]*"( P-value )"))


## Volcano plot of P-value vs Survival Difference
library(ggplot2)
library(ggrepel)
plot_volcano_survival = function(df, days_cutoff = 2, P_cutoff = 0.05, OS_type = T) {
  # df = dataframe containing plot data
  # days_cutoff = threshold of biological significance
  # P_cutoff = threshold of statistical significance
  # OS_type = T is "OS" and F is "PFS"
  if (OS_type) {
    df$log.pvalue = -log10(df$OS_pval)
    df$difference = df$median_OS_diff
  } else {
    df$log.pvalue = -log10(df$PFS_pval)
    df$difference = df$median_PFS_diff
  }
  
  df$color = "black"
  df$color[df$Gene %in% c("HSPA9", "HSPA5", "HSPA4", "HSPA8", "HSPA13")] = "blue"
  df$color[grepl("PSMA", df$Gene)] = "red"
  df$color[grepl("PSMB", df$Gene)] = "darkgreen"
  df$color[grepl("PSMC|PSMD", df$Gene)] = "purple"
  df$trans = 0.1   # point transparency
  df$trans[df$color != "black"] = 0.7
  df$color = factor(df$color)
  
  fig = ggplot(df) +
    geom_point(aes(difference, log.pvalue, colour = color), alpha = df$trans) +
    geom_point(data = subset(df, color != 'black'),
               aes(x = difference, 
                   y = log.pvalue, 
                   colour = color),
               alpha = subset(df, color != "black")$trans) +
    scale_colour_manual(values = levels(df$color)) + 
    geom_hline(yintercept = -log10(P_cutoff), linetype = "dashed", color = "gray50") +
    geom_vline(xintercept = c(days_cutoff, -days_cutoff), linetype = "dashed", color = "gray50") +
    scale_x_continuous(
      name = expression("Difference in PFS Survival in Days (Top 20% - Bottom 20%)")) +
    scale_y_continuous(
      name = expression("- log"[10]*"( P-value )")) +
    theme_minimal() +
    theme(legend.position = "none") + 
    geom_text_repel(
      data = subset(df, df$Gene %in% c("HSPA9", "HSPA5", "HSPA4", "HSPA8", "HSPA13", 
                                       #"PSMC2", "PSMD2",
                                       grep("PSMC|PSMD", df$Gene, value = T),
                                       grep("PSMA", df$Gene, value = T),
                                       grep("PSMB", df$Gene, value = T))),
      aes(difference, log.pvalue, label = Gene, colour = color),
      size = 3)
  
  fig
}
plot_volcano_survival(x, days_cutoff = 100, P_cutoff = 0.05, OS_type = F)  # 6x4


