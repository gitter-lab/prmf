#!/usr/bin/env Rscript
library(argparse)
library(CoGAPS)

main = function() {
  parser = ArgumentParser(description="Run CoGAPS")
  parser$add_argument('--data', required=TRUE)
  parser$add_argument('--k-latent', required=TRUE)
  parser$add_argument('--outdir', required=TRUE)
  args = parser$parse_args()

  data = as.matrix(read.csv(args$data, row.names=1))
  params = new("CogapsParams")
  params = setParam(params, "nPatterns", args$k_latent)

  # data = AP
  results = CoGAPS(data, params)
  write.csv(results$Amean, file.path(args$outdir, 'A.csv'), row.names=FALSE, col.names=FALSE)
  write.csv(results$Pmean, file.path(args$outdir, 'P.csv'), row.names=FALSE, col.names=FALSE)
}

main()
