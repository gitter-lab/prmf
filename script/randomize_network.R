#!/usr/bin/env Rscript
library(BiRewire)
library(igraph)
library(argparse)

main = function() {
  parser = ArgumentParser(description='Randomize edges in a network in a degree-controlled manner')
  parser$add_argument('--in-network', help='Network to randomize in graphml format', required=T)
  parser$add_argument('--out-network', help='Randomized network', required=T)
  args = parser$parse_args()

  G = read_graph(args$in_network, format='graphml')
  H = birewire.rewire.undirected(G)
  write_graph(H, args$out_network, format='graphml')
}

main()
