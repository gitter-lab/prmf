#!/usr/bin/env python
import sys, argparse
import os, os.path
import networkx as nx

def main():
  parser = argparse.ArgumentParser(description="""
For each graphml file in indir, write a nodelist in outdir with the same basename
""")
  parser.add_argument("--indir", type=str, required=True)
  parser.add_argument("--outdir", type=str, required=True)
  args = parser.parse_args()

  for bn_ext in os.listdir(args.indir):
    bn, ext = os.path.splitext(bn_ext)
    if ext == ".graphml":
      ifp = os.path.join(args.indir, bn_ext)
      ofp = os.path.join(args.outdir, bn + ".txt")
      G = nx.read_graphml(ifp)
      with open(ofp, 'w') as ofh:
        ofh.write("\n".join(sorted(G.nodes())))

if __name__ == "__main__":
  main()
