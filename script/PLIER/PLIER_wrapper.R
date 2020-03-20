#!/usr/bin/env Rscript
library(PLIER)
library(argparse)

# Embed vec into <dimen>-dimensional space
embed_vec = function(dimen, vec, ind_map) {
  rv = array(data = 0, dim=dimen)
  for (i in 1:length(vec)) {
    j = ind_map[i]
    rv[j] = vec[i]
  }
  rv
}

# Embed columns of <arr> into <dimen>-dimensional space
# <ind_map> associates a row in the old vector space to a row in the new one
embed_array = function(dimen, arr, ind_map) {
  dim_v = dim(arr)
  n_row = dim_v[[1]]
  n_col = dim_v[[2]]
  rv = array(data=NA, dim=c(dimen, n_col))
  for (j in 1:n_col) {
    rv[,j] = embed_vec(dimen, arr[,j], ind_map)
  }
  rv
}

embed_array_named = function(arr, row_names) {
  dim_v = dim(arr)
  n_col = dim_v[[2]]
  rv = array(data=0, dim=c(length(row_names), n_col))
  row.names(rv) = row_names
  rv[row.names(arr),] = arr
  rv
}

# Return a integer vector which associates 
map_indexes = function(all_names, some_names) {
  rv = 
  rv
}

# Create binary 1-column matrix from a character vector
to_binary_array = function(vec) {
  rv = matrix(rep(1, length(vec)), nrow=length(vec))
  rownames(rv) = vec
  rv
}

# Split file extension off filepath or filename
splitext = function(text) {
  rv = NA
  pattern = "([^.]*)([.][^.]*)$"
  match_data = regmatches(text, regexec(pattern, text))
  if(length(match_data[[1]]) != 1) {
    rv = list(match_data[[1]][2], match_data[[1]][3])
  }
  rv
}

main2 = function(args) {
  # assume data has row and column names
  # assume data is such that there are more features than observations and that features belong on rows
  data = read.csv(args$data, row.names=1)
  data_dim = dim(data)
  if(data_dim[1] < data_dim[2]) {
    data = t(data)
    data_dim = dim(data)
  }
  # If Z-score normalization, error is thrown:
  # Error in svd(Y, nu = k, nv = k) : infinite or missing values in 'x'
  data = rowNorm(data)

  # PLIER initialization requires that all genes have some measurement
  # subset the data to meet this requirement
  data_rowSums = rowSums(data)
  inds = which(data_rowSums != 0)
  warning(paste0("Subsetting data matrix to remove rows without measurements: ", data_dim[1], " -> ", length(inds)))
  data = data[inds,]
  data_dim = dim(data)
  
  prior = NULL
  pathway_names = list()
  if(!is.null(args$pathways_file)) {
    library(igraph)
    pathways = list()
    con = file(args$pathways_file, "r")
    i = 0
    while ( TRUE ) {
      line = readLines(con, n = 1)
      if ( length(line) == 0 ) {
        break
      }
      i = i + 1
      pathway_graph = read.graph(line, format='graphml')
      pathway_matrix = matrix(rep(1, vcount(pathway_graph)), nrow=vcount(pathway_graph))
      rownames(pathway_matrix) = vertex_attr(pathway_graph, args$node_attribute, index = V(pathway_graph))
      bn_ext = basename(line)
      splitext_rv = splitext(bn_ext)
      bn = splitext_rv[[1]]
      ext = splitext_rv[[2]]
      colnames(pathway_matrix) = bn
      pathways[[i]] = pathway_matrix
      pathway_names[[i]] = bn
    }
    close(con)
    prior = do.call(combinePaths, pathways)
    colnames(prior) = pathway_names
  } else {
    tbl = read.table(args$pathways_csv, sep=',')
    row.names = rownames(tbl)
    prior = as.matrix(tbl)
    rownames(prior) = row.names
    pathway_names = colnames(prior)
  }
  
  # TODO add colnames to be able to pass computeAUC=T
  rownames_inter = commonRows(data, prior) # "common genes"
  rownames_missing = setdiff(row.names(data), rownames_inter)

  if(args$extra_prior) {
    warning(paste0("Adding extra prior to include ", length(rownames_missing), " rows in data not found in any other prior"))
    extra_prior = to_binary_array(rownames_missing)
    prior = combinePaths(prior, extra_prior)
    pathway_names[[i+1]] = 'extra'
  } else {
    warning(paste0("Subsetting data matrix to only include rows shared with a prior: ", data_dim[1], " -> ", length(rownames_inter)))
  }
  
  # from source, return vaue is:
  # list(residual=(Y-Z%*%B), B=B, Z=Z, U=U, C=C, numActPath=length(ii), L1=L1, L2=L2, L3=L3, heldOutGenes=heldOutGenes)
  # TODO numActPat, heldOutGenes?
  # TODO other things added to namespace if computeAUC
  # TODO remove minGenes
  write('colnames', stdout())
  write(colnames(prior), stdout())
  # NOTE if PLIER fails because B is ill-conditioned, increasing L1 and/or L2 may eliminate the problem due to how the update rules for B are defined
  # NOTE if args$L1 or args$L2 are null, PLIER will generate defaults
  plierResult = PLIER(as.matrix(data[rownames_inter,]), prior[rownames_inter,], k=args$k_latent, trace=T, seed=args$seed, L1=args$L1, L2=args$L2)
  write.csv(plierResult$residual, file.path(args$outdir, 'residual.csv'))
  write.csv(plierResult$B, file.path(args$outdir, 'B.csv'))
  write.table(plierResult$Z, file.path(args$outdir, 'Z.csv'), col.names=F, row.names=F, sep=",")

  # add pathway names to U as well (just as for C := prior)
  # Z - CU
  # Z is gene x latent
  # C is gene x pathway
  # U is pathway x latent
  row.names(plierResult$U) = pathway_names
  write.csv(plierResult$U, file.path(args$outdir, 'U.csv'))

  write.csv(plierResult$C, file.path(args$outdir, 'C.csv'))
  con = file(file.path(args$outdir, 'opt.txt'), 'w')
  writeLines(c(
    sprintf('data_residual = %f', norm(plierResult$residual, type='F')),
    sprintf('prior_residual = %f', norm(plierResult$Z - plierResult$C %*% plierResult$U)),
    sprintf('L1 = %f', plierResult$L1), 
    sprintf('L2 = %f', plierResult$L2),
    sprintf('L3 = %f', plierResult$L3)
  ), con)
  close(con)
}

# R's argparse doesnt implement add_mutually_exclusive_group(required=T)
is_mutually_exclusive_required_group = function(env, members) {
  rv = F
  env = as.environment(env)
  n_not_null = 0
  for (i in 1:length(members)) {
    member_v = get(members[[i]], envir=env)
    if (!is.null(member_v)) {
      n_not_null = n_not_null + 1
    }
  }
  if (n_not_null == 1) {
    rv = T
  }
  return(rv)
}

main = function() {
  parser = ArgumentParser(description='Construct binary matrices from pathway files and run PLIER')
  parser$add_argument('--data', required=T)
  parser$add_argument('--pathways-csv')
  parser$add_argument('--pathways-file')
  parser$add_argument('--outdir', required=T)
  parser$add_argument('--k-latent', default=6, type='integer')
  parser$add_argument('--L1', help='Regularization parameter for reconstructing Z (the PLIER-generated gene by latent matrix) with the product of C (the user-provided prior matrix) and U (the PLIER-generated prior by latent matrix). If not provided, use PLIER default', type='double')
  parser$add_argument('--L2', help='Regularization parameter for Frobenius norm on B (the PLIER-generated latent by sample matrix): the shrinkage parameter on B. If not provided, use PLIER default', type='double')
  parser$add_argument('--seed', default=1, type='integer')
  parser$add_argument('--node-attribute', default='id', help="graphml node attribute to use for gene identifiers")
  parser$add_argument('--extra-prior', action='store_true')
  args = parser$parse_args()
  if (!is_mutually_exclusive_required_group(args, c('pathways_csv', 'pathways_file'))) {
    writeLines('Exactly one of --pathways-csv or --pathways-file is required', con=stderr())
    q(status=2)
  }
  main2(args)
}

main()
