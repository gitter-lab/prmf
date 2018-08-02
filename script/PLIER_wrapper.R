#!/usr/bin/env Rscript
library(PLIER)
library(igraph)
library(argparse)

main = function() {
  parser = ArgumentParser(description="
Construct binary matrices from pathway files and run PLIER

PLIER_wrapper.R [-h] data.csv pathways_file.txt
")
  parser$add_argument('--data', required=T)
  parser$add_argument('--nodelist', required=T)
  parser$add_argument('--pathways-file', required=T)
  parser$add_argument('--outdir', required=T)
  parser$add_argument('--k-latent', default=6, type='integer')
  args = parser$parse_args()

  # attach row names to data
  data = read.csv(args$data)
  row.names(data) = readLines(args$nodelist)

  pathways = list()
  con = file(args$pathways_file, "r")
  i = 0
  while ( TRUE ) {
    line = readLines(con, n = 1)
    if ( length(line) == 0 ) {
      break
    }
    i += 1
    pathway_graph = read.graph(line, format='graphml')
    pathway_matrix = matrix(rep(1, vcount(pathway_graph)), nrow=vcount(pathway_graph))
    rownames(pathway_matrix) = vertex_attr(pathway_graph, 'id', index = V(pathway_graph))
    pathways[[i]] = pathway_matrix
  }
  close(con)

  prior = combinePaths(pathways...)
  all_rownames = commonRows(data, prior) # "common genes"

  # from source, return vaue is:
  # list(residual=(Y-Z%*%B), B=B, Z=Z, U=U, C=C, numActPath=length(ii), L1=L1, L2=L2, L3=L3, heldOutGenes=heldOutGenes)
  # TODO numActPat, heldOutGenes?
  plierResult = PLIER(data[all_rownames,], prior[all_rownames,], k = args$k_latent, trace = T)
  write.csv(plierResult$residual, file.path(args$outdir, 'residual.csv'))
  write.csv(plierResult$B, file.path(args$outdir, 'B.csv'))
  write.csv(plierResult$Z, file.path(args$outdir, 'Z.csv'))
  write.csv(plierResult$U, file.path(args$outdir, 'U.csv'))
  write.csv(plierResult$C, file.path(args$outdir, 'C.csv'))
  con = file(file.path(args$outdir, 'opt.txt'), 'w')
  writeLines(c(
    sprintf('L1 = ', plierResult$L1), 
    sprintf('L2 = ', plierResult$L2),
    sprintf('L3 = ', plierResult$L3)
  ), con)
}

main()
