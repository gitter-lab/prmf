#!/usr/bin/env Rscript
library(argparse)
library(CoGAPS)

main = function() {
  parser = ArgumentParser(description="Run CoGAPS")
  parser$add_argument('-d', '--data', required=TRUE, help="Data matrix with shape n_samples x n_features")
  parser$add_argument('-k', '--k-latent', required=FALSE, type='integer', default=7, help="CoGAPS nPatterns parameter")
  parser$add_argument('-o', '--outdir', required=TRUE, help="Location to write results")
  parser$add_argument('-t', '--transpose-data', required=FALSE, action='store_true', default=FALSE, help="Provide this flag if the data file is n_features x n_samples")
  parser$add_argument('-y', '--data-type', required=FALSE, default='csv', help="Provide --data-type RDS if <--data> is an RDS file rather than a CSV")
  args = parser$parse_args()

  # Load data
  write('[STATUS] Loading data...', stdout())
  data = NULL
  if (args$data_type == 'csv') {
    data = as.matrix(read.csv(args$data, row.names=1))
  } else if (args$data_type == 'RDS') {
    data = readRDS(args$data)
  } else {
    write(sprintf('Unrecognized --data-type = %s', args$data_type), stderr())
    quit('save'='no', status=2)
  }
  dim_rv = dim(data)
  write(sprintf('[STATUS] Done loading data with shape = (%d, %d)', dim_rv[1], dim_rv[2]), stdout())

  # Check for non-numeric row names in the first column
  if (class(data[2,1]) == 'character') {
    row.names(data) = data[,1]
    data = data[,2:dim_rv[2]]
  }

  # Prepare CoGAPS parameters
  params = new("CogapsParams")
  params = setParam(params, "nPatterns", args$k_latent)

  # Run CoGAPS
  # this script expects data to be samples x genes
  # CoGAPS requires data to be genes x samples
  # sampleFactors is samples x latent
  # featureLoadings is feature x latent
  # data \approx sampleFactors \cdot featureLoadings
  results = CoGAPS(data, params, transposeData=!args$transpose_data, outputFrequency=10)
  write.csv(results@sampleFactors, file.path(args$outdir, 'sample_by_latent.csv'))
  write.csv(results@featureLoadings, file.path(args$outdir, 'feature_by_latent.csv'))
  return(results)
}
results = main()
# there is some metadata that could be written, interrogate with save.image() if desired
