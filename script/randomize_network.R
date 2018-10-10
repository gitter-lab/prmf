#!/usr/bin/env Rscript
load_success = TRUE
load_success = load_success && library(BiRewire, quietly=TRUE, logical.return=TRUE)
load_success = load_success && library(igraph, quietly=TRUE, logical.return=TRUE)
load_success = load_success && library(argparse, quietly=TRUE, logical.return=TRUE)
if(!load_success) {
  write("Unable to load dependencies", stderr())
  quit(status=2)
}

# Split file extension off filepath or filename
# TODO copied from PLIER_wrapper.R
splitext = function(text) {
  rv = NA
  pattern = "([^.]*)([.][^.]*)$"
  match_data = regmatches(text, regexec(pattern, text))
  if(length(match_data[[1]]) != 1) {
    rv = list(match_data[[1]][2], match_data[[1]][3])
  }
  rv
}

main = function() {
  parser = ArgumentParser(description='Randomize edges in a network in a degree-controlled manner')
  parser$add_argument('--indir', help='Network to randomize in graphml format', required=T)
  parser$add_argument('--outdir', help='Randomized network', required=T)
  args = parser$parse_args()

  ll = list.files(args$indir)
  for (ifn in ll) {
    is_graphml = FALSE
    ifp = file.path(args$indir, ifn)
    splitext_rv = splitext(basename(ifn))
    if (!is.na(splitext_rv)) {
      bn = splitext_rv[1]
      ext = splitext_rv[2]
      if (ext == '.graphml') {
        is_graphml = TRUE
      }
    }

    if (is_graphml) {
      ofp = file.path(args$outdir, ifn)

      G = read_graph(ifp, format='graphml')
      H = birewire.rewire.undirected(G)
      write_graph(H, ofp, format='graphml')
    }
  }
}

main()
