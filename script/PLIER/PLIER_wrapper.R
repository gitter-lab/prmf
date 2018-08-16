#!/usr/bin/env Rscript
library(PLIER)
library(igraph)
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

splitext = function(text) {
  rv = NA
  pattern = "([^.]*)([.][^.]*)$"
  match_data = regmatches(text, regexec(pattern, text))
  if(length(match_data) != 1) {
    rv = list(match_data[[1]][2], match_data[[1]][3])
  }
  rv
}

main2 = function(args) {
  # attach row names to data
  # assume data is such that there are more features than observations and that features belong on rows
  data = read.csv(args$data, header=F)
  dim_data = dim(data)
  if(dim_data[1] < dim_data[2]) {
    data = t(data)
  }
  # If Z-score normalization, error is thrown:
  # Error in svd(Y, nu = k, nv = k) : infinite or missing values in 'x'
  data = rowNorm(data)
  nodelist = readLines(args$nodelist)
  row.names(data) = nodelist

  # PLIER initialization requires that all genes have some measurement
  # subset the data to meet this requirement
  data_rowSums = rowSums(data)
  inds = which(data_rowSums != 0)
  data = data[inds,]
  
  pathways = list()
  pathway_names = list()
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
    rownames(pathway_matrix) = vertex_attr(pathway_graph, 'id', index = V(pathway_graph))
    bn_ext = basename(line)
    splitext_rv = splitext(bn_ext)
    bn = splitext_rv[1]
    ext = splitext_rv[2]
    colnames(pathway_matrix) = bn
    pathways[[i]] = pathway_matrix
    pathway_names[[i]] = bn
  }
  close(con)
  
  # TODO add colnames to be able to pass computeAUC=T
  prior = do.call(combinePaths, pathways)
  all_rownames = commonRows(data, prior) # "common genes"
  
  # from source, return vaue is:
  # list(residual=(Y-Z%*%B), B=B, Z=Z, U=U, C=C, numActPath=length(ii), L1=L1, L2=L2, L3=L3, heldOutGenes=heldOutGenes)
  # TODO numActPat, heldOutGenes?
  # TODO other things added to namespace if computeAUC
  plierResult = PLIER(data[all_rownames,], prior[all_rownames,], k=args$k_latent, trace=T, computeAUC=F, seed=args$seed)
  write.csv(plierResult$residual, file.path(args$outdir, 'residual.csv'))
  write.csv(plierResult$B, file.path(args$outdir, 'B.csv'))

  # embed Z in #{nodelist}-dimensional space (use association above)
  Z_embed = embed_array(length(nodelist), plierResult$Z, inds)
  write.table(Z_embed, file.path(args$outdir, 'Z.csv'), col.names=F, row.names=F, sep=",")

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

main = function() {
  parser = ArgumentParser(description='Construct binary matrices from pathway files and run PLIER')
  parser$add_argument('--data', required=T)
  parser$add_argument('--nodelist', required=T)
  parser$add_argument('--pathways-file', required=T)
  parser$add_argument('--outdir', required=T)
  parser$add_argument('--k-latent', default=6, type='integer')
  parser$add_argument('--seed', default=1, type='integer')
  args = parser$parse_args()
  main2(args)
}

main()
