#!/usr/bin/env Rscript
library(argparse)
library(tidyr)

parser = ArgumentParser(description="Apply MultiPLIER-style transfer learning")
parser$add_argument('--gene-by-sample-tidy', '-y', help='Data for Y in the form of a tidy dataset', required=TRUE)
parser$add_argument('--tidy-row', '-r', required=TRUE, help="Tidy variable to use as rows in matrix")
parser$add_argument('--tidy-col', '-c', required=TRUE, help="Tidy variable to use as columns in matrix")
parser$add_argument('--tidy-val', '-v', required=TRUE, help="Tidy variable to use as values in matrix")
parser$add_argument('--tidy-meta', '-m', required=TRUE, help="Meta data for each row/sample e.g. tumorType")
parser$add_argument('--gene-by-latent', '-z', help='A learned gene-by-latent matrix from a prior run of a matrix factorization method e.g. PRMF, CoGAPS, PLIER, etc.', required=TRUE)
parser$add_argument('--outdir', '-o', help='Directory to write results to', required=TRUE)
args = parser$parse_args()

# PLIER: Y = ZB
# Y is gene x sample
# Z is gene x latent
# B is latent x sample
y_tidy = read.csv(args$gene_by_sample_tidy, sep='\t')
y_df_t = y_tidy %>% 
  dplyr::group_by(!!as.name(args$tidy_row), !!as.name(args$tidy_col)) %>% 
  dplyr::summarise(mean = mean(!!as.name(args$tidy_val))) %>% 
  spread(!!as.name(args$tidy_col), mean)
sample_names = levels(y_df_t[[args$tidy_row]])

# write sample metadata
sample_meta = y_tidy %>% 
  dplyr::select(!!as.name(args$tidy_row), !!as.name(args$tidy_meta)) %>% 
  dplyr::group_by(!!as.name(args$tidy_row)) %>% 
  dplyr::summarise(first = dplyr::first(!!as.name(args$tidy_meta))) %>% 
  dplyr::mutate(!!as.name(args$tidy_meta) := first) %>%
  dplyr::select(!!as.name(args$tidy_row), !!as.name(args$tidy_meta))
write.csv(sample_meta, file.path(args$outdir, 'sample_meta.csv'), row.names=FALSE)

# Z has row names for the genes
z_matrix = read.csv(args$gene_by_latent, row.names=1)
z_row_names = row.names(z_matrix)
z_matrix = data.matrix(z_matrix)

# Z may be defined on ENSG or HGNC
# map ENSG to HGNC if needed
# map identifiers using HGNC custom download
if (substr(row.names(z_matrix)[1], 1, 4) == "ENSG") {
  url = 'https://www.genenames.org/cgi-bin/download/custom?col=gd_hgnc_id&col=gd_app_sym&col=gd_app_name&col=gd_status&col=gd_prev_sym&col=gd_aliases&col=gd_pub_chrom_map&col=gd_pub_acc_ids&col=gd_pub_refseq_ids&col=gd_pub_ensembl_id&status=Approved&status=Entry%20Withdrawn&hgnc_dbtag=on&order_by=gd_app_sym_sort&format=text&submit=submit'
  ensg_hgnc_map_fp = file.path(args$outdir, 'ensg_hgnc_map.tsv')
  if (!file.exists(ensg_hgnc_map_fp)) {
    download.file(url, ensg_hgnc_map_fp)
  }
  ensg_hgnc_df = read.csv(ensg_hgnc_map_fp, sep='\t')
  ensg_hgnc_map = new.env()
  build_map = function(x) {
    if (x[['Ensembl.gene.ID']] != "") {
      ensg_hgnc_map[[x[['Ensembl.gene.ID']]]] = x[['Approved.symbol']]
    }
  }
  apply_rv = apply(ensg_hgnc_df, 1, build_map)
  new_row_names = sapply(z_row_names, function(ensg) {
    rv = ensg
    hgnc = ensg_hgnc_map[[ensg]]
    if(!is.null(hgnc)) {
      rv = hgnc
    }
    return(rv)
  })
  new_row_names = as.vector(new_row_names)
  row.names(z_matrix) = new_row_names
}

# Report initial numbers of genes before intersection
y_matrix = t(data.matrix(y_df_t))
y_matrix[is.na(y_matrix)] = 0
write(paste0("Number of genes in input sample x gene matrix Y: ", dim(y_matrix)[1]), stderr())
write(paste0("Number of genes in input gene x latent matrix Z: ", dim(z_matrix)[1]), stderr())

# Only select genes that the transfer learning model was trained on
common_genes = intersect(row.names(z_matrix), colnames(y_df_t))
write(paste0("Number of genes common to input datasets after mapping: ", length(common_genes)), stderr())
if (length(common_genes) == 0) {
  print("No common genes, quitting")
  quit(status=1)
}
y_matrix = y_matrix[common_genes,]
z_matrix = z_matrix[common_genes,]
  
# L2 defaults to smallest singular value in SVD according to PLIER.
# it's probably the case that this parameter does not affect differential expression analysis using B.
# however, this is a poor default (especially in low sample settings): the smallest singular 
# value will most likely be a small value and make the L2 regularizer irrelevant.
# it is a tradeoff between the reconstruction error and the magnitude of values in B.
# these quantities are proportional to the size of Y and B respectively (by virtue of the Frobenius norm).
# therefore, the L2 parameter should also be on that scale to be appreciable.
# for example, the smallest singular value from the recount2 data is 1.6, and in some experiments the
# effect of the L2 parameter is negligible at least up to a setting of 1100 (the ratio of the number 
# of entries in Y to the number of entries in B): the change in the Frobenius norm of B is 
# less than a half a percent. Settings should be a multiplicative factor of this ratio to take effect.
# A 2% decrease in the norm of B was observed with a factor of 10, while a 17% decrease was observed with
# a factor of 100.
m_genes = dim(y_matrix)[1]
n_samples = dim(y_matrix)[2]
k_latent = dim(z_matrix)[2]
ratio = (m_genes * n_samples) / (k_latent * n_samples)
L2 = ratio * 100

# Apply PLIER model to solve for B
# https://github.com/wgmao/PLIER/blob/a2d4a2aa343f9ed4b9b945c04326bebd31533d4d/R/Allfuncs.R#L465
# B = (Z^T Z + \lambda * I)^{-1} * Z^T Y
# this is the formula for a pseudoinverse of Z applied to Y but where there is a L2 regularizer on B
b_matrix = solve(t(z_matrix)%*%z_matrix+L2*diag(k_latent))%*%t(z_matrix)%*%y_matrix
colnames(b_matrix) = sample_names
write.csv(t(b_matrix), file.path(args$outdir, 'sample_by_latent_transfer.csv'))
