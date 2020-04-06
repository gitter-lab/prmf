#!/usr/bin/env Rscript
requireNamespace('argparse', quietly=TRUE)
library(ggplot2)
library(dplyr)

parser = argparse::ArgumentParser(description="Run Kolmogorov-Smirinov Test for difference in distributions of latent variable loading across tumor types")
parser$add_argument("--b-matrix", "-b", required=TRUE)
parser$add_argument("--sample-meta", "-s", required=TRUE)
parser$add_argument("--outdir", "-o", required=TRUE)
args = parser$parse_args()

# NF acronyms
# TODO what if some other tumorType is present in --sample-meta
abbrev_df = data.frame(
  tumorType=c("Cutaneous Neurofibroma", "High Grade Glioma", "Low Grade Glioma", "Malignant Peripheral Nerve Sheath Tumor", "Meningioma", "Neurofibroma", "Plexiform Neurofibroma", "Schwannoma", ""),
  tumorTypeAbbrev=c('cNF', 'HGG', 'LGG', 'MPNST', 'MG', 'NF', 'pNF', 'SW', "")
)

n_pairwise_plots = 1
n_per_class_plots = 1

b_matrix = read.csv(args$b_matrix)
# set colnames for df joins
b_matrix_cols = colnames(b_matrix)
b_matrix_cols[1] = 'specimenID'
colnames(b_matrix) = b_matrix_cols

n_samples = dim(b_matrix)[1]
k_latent = dim(b_matrix)[2] - 1
sample_meta = read.csv(args$sample_meta)

# count number of samples of each tumor type
tumor_count_df = sample_meta %>%
  group_by(tumorType) %>%
  summarise(count=dplyr::n())
abbrev_and_count_df = merge(x=abbrev_df, y=tumor_count_df, by='tumorType')

# Conduct KS-test for each pair of tumor types for each LV
tumor_types = levels(sample_meta$tumorType)
n_tests = length(tumor_types) * k_latent
ks_summary_tumor = c()
ks_summary_lv = c()
ks_summary_p = c()
for (tumor_type in tumor_types) {
  for (k in 2:(k_latent+1)) { # +1 because first column of b_matrix is the sample name column
    x = b_matrix[sample_meta$tumorType == tumor_type, k]
    y = b_matrix[sample_meta$tumorType != tumor_type, k]
    ks_rv = ks.test(x, y)
    p_value = p.adjust(ks_rv$p.value, method='bonferroni', n=n_tests)

    ks_summary_tumor = c(ks_summary_tumor, tumor_type)
    ks_summary_lv = c(ks_summary_lv, k - 1) # - 1 for the real LV number
    ks_summary_p = c(ks_summary_p, p_value)
  }
}
ks_summary = data.frame(tumorType=ks_summary_tumor, LV=ks_summary_lv, pValue=ks_summary_p)

# Plot the best LV which differentiates one tumor type from the rest for each tumor type
ks_summary$rowNum = 1:dim(ks_summary)[1]
ks_summary_group_tumor = ks_summary %>%
  group_by(tumorType) %>%
  slice(which.min(pValue))
ks_summary_group_tumor = as.data.frame(ks_summary_group_tumor)
for (i in 1:dim(ks_summary_group_tumor)[1]) {
  k = ks_summary_group_tumor[i,]$LV
  k_label = colnames(b_matrix)[k+1]
  tumor_type = ks_summary_group_tumor[i,]$tumorType
  if (ks_summary_group_tumor[i,'pValue'] == 0) {
    p_value_str = sprintf("KS-test p < %1.3e", .Machine$double.eps)
  } else {
    p_value_str = sprintf('KS-test p = %1.3e', ks_summary_group_tumor[i,'pValue'])
  }

  tumor_type_data = abbrev_and_count_df[abbrev_and_count_df$tumorType == tumor_type,]
  tumor_type_abbrev = tumor_type_data$tumorTypeAbbrev

  val_x = b_matrix[sample_meta$tumorType == tumor_type, k_label]
  label_x = paste0(tumor_type_abbrev, sprintf('\nN = %d', length(val_x)))
  label_x_r = replicate(length(val_x), label_x)

  val_y = b_matrix[sample_meta$tumorType != tumor_type, k_label]
  label_y = paste0('not ', tumor_type_abbrev, sprintf('\nN = %d', length(val_y)))
  label_y_r = replicate(length(val_y), label_y)

  plot_df = data.frame(val=c(val_x, val_y), label=c(label_x_r, label_y_r))
  plot = ggplot(plot_df, aes(x=label, y=val)) +
    geom_boxplot() +
    ggtitle(paste0("Differential Loading\n", p_value_str)) + 
    xlab("Tumor type") +
    ylab(paste0("LV", k, ' loading')) + theme(
    plot.title = element_text(size=24, hjust=0.5),
    axis.title.x = element_text(size=20),
    axis.title.y = element_text(size=20),
    axis.text = element_text(size=16)
    )
  ggsave(file.path(args$outdir, paste0('plot_', tumor_type, '.png')), plot=plot)
}

# Summarize pairwise tests for each LV into a single metric
ks_summary_group_lv = ks_summary %>% 
  group_by(LV) %>% 
  summarise(mean=mean(pValue), median=median(pValue), min=min(pValue), max=max(pValue))
ks_summary_group_lv_order = order(ks_summary_group_lv$median)

# Plot the most differentially expressed LV according to the median KS p-value
for (i in 1:n_per_class_plots) {
  row = ks_summary_group_lv_order[i]
  k = ks_summary_group_lv[row,]$LV
  k_label = colnames(b_matrix)[k+1]
  df_merge = merge(x=b_matrix, y=sample_meta, by='specimenID')
  df_merge = merge(x=df_merge, y=abbrev_df, by='tumorType')
  df_merge = merge(x=df_merge, y=tumor_count_df, by='tumorType')
  df_merge = df_merge %>%
    dplyr::mutate(tumorTypeCount=sprintf('%s (%d)', tumorTypeAbbrev, count))
  plot = ggplot(df_merge, aes(x=tumorTypeCount, y=!!as.name(k_label))) +
    geom_boxplot() +
    ggtitle(paste0("Differential Loading per Class")) + 
    xlab("Tumor type") +
    ylab(paste0("LV", k, ' loading')) + theme(
    plot.title = element_text(size=24, hjust=0.5),
    axis.text.x = element_text(size=16, angle=90, hjust=1, vjust=0.5),
    axis.title.y = element_text(size=20),
    axis.title.x = element_text(size=20),
    axis.text.y = element_text(size=16)
    )
  ggsave(file.path(args$outdir, paste0('plot_per_class', i, '.png')), plot=plot)
}
