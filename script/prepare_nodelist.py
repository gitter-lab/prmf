#!/usr/bin/env python
import sys, argparse
import os, os.path
import networkx as nx
import factorlib as fl
from factorlib import string_db as sdb

def main():
  parser = argparse.ArgumentParser()
  parser.add_argument("--stringdb", type=argparse.FileType('r'), required=True)
  parser.add_argument("--graphml-dir", type=str, required=True)
  parser.add_argument("--out-nodelist", type=argparse.FileType('w'), required=True)
  parser.add_argument("--out-graph", type=argparse.FileType('wb'), required=True)
  args = parser.parse_args()

  # TODO use fl.prepare_nodelist for consistency
  Gs = fl.parse_pathways(args.graphml_dir)
  G_ppi = sdb.parse_string_fh(args.stringdb)
  Gs = [G_ppi] + Gs

  G_union = fl.weighted_union(Gs)
  nodes = sorted(G_union.nodes())
  for node in nodes:
    args.out_nodelist.write(node + "\n")

  nx.write_graphml(G_union, args.out_graph)

if __name__ == "__main__":
  main()
