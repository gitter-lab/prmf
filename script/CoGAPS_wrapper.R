#!/usr/bin/env Rscript
library(argparse)
library(CoGAPS)

main = function() {
  parser = ArgumentParser(description="Run CoGAPS")
  parser$add_argument('--data')
  parser$add_argument('--k-latent')
  parser$add_argument('--outdir')
  args = parser$parse_args()

  D = as.matrix(read.csv(args$data, header=FALSE))
  n_obs = dim(D)[1]
  n_features = dim(D)[2]

  # generate stdev measurements by computing the standard deviation for each observation
  # and using that as the stdev for each feature measurement in that observation
  row_stdev = apply(D, 1, sd)
  S = as.matrix(do.call(cbind, replicate(n_features, row_stdev, simplify=FALSE)))

  # D = AP
  results = CoGAPS(D, S, nFactor = as.numeric(args$k_latent))
  write.csv(results$Amean, file.path(args$outdir, 'A.csv'), row.names=FALSE, col.names=FALSE)
  write.csv(results$Pmean, file.path(args$outdir, 'P.csv'), row.names=FALSE, col.names=FALSE)
}

main()
