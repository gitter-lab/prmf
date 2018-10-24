#!/usr/bin/env python
import sys, argparse
import os, os.path
import networkx as nx
import factorlib as fl
from factorlib import string_db as sdb

def main():
  parser = argparse.ArgumentParser(description="""
Construct a nodelist containing all nodes from STRING and all nodes from all networks in <graphml-dir>.
""")
  parser.add_argument("--stringdb", type=argparse.FileType('r'), required=True)
  parser.add_argument("--graphml-dir", type=str)
  parser.add_argument("--graphmls", type=str, nargs='+')
  parser.add_argument("--node-attribute", type=str, help='Node attribute which contains protein identifier; if not provided, the node identifier is used as a protein identifier')
  parser.add_argument("--out-nodelist", type=argparse.FileType('w'), required=True)
  parser.add_argument("--out-graph", type=argparse.FileType('wb'), required=True)
  args = parser.parse_args()

  if args.graphml_dir is None and args.graphmls is None:
    sys.stderr.write("Exactly one of --graphml-dir or --graphmls is required.\n")
    sys.exit(21)
  if args.graphml_dir is not None and args.graphmls is not None:
    sys.stderr.write("Exactly one of --graphml-dir or --graphmls is required.\n")
    sys.exit(22)

  # TODO use fl.prepare_nodelist for consistency
  Gs = []
  if args.graphmls is not None:
    Gs = fl.parse_pathways(args.graphmls)
  if args.graphml_dir is not None:
    Gs = fl.parse_pathways_dir(args.graphml_dir)

  # relabel nodes if needed
  Gs = list(map(lambda G: fl.relabel_nodes(G, args.node_attribute), Gs))

  G_ppi = sdb.parse_string_fh(args.stringdb)

  # check identifiers, error if there is no overlap with --stringdb
  node_ids = set()
  for G in Gs:
    for node in G.nodes():
      node_ids.add(node)
  all_missing = True
  for node in node_ids:
    if node in G_ppi:
      all_missing = False
      break
  if(all_missing):
    sys.stderr.write("No node identifiers in the provided pathways overlap with the background protein-protein interaction network. Exiting.\n")
    sys.exit(23)

  Gs = [G_ppi] + list(Gs)

  G_union = fl.weighted_union(Gs)
  nodes = sorted(G_union.nodes())
  for node in nodes:
    args.out_nodelist.write(node + "\n")

  nx.write_graphml(G_union, args.out_graph)

if __name__ == "__main__":
  main()
