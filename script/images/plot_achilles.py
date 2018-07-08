#!/usr/bin/env python
import argparse 
import os, os.path
import matplotlib.pyplot as plt
import numpy as np
import factorlib as fl

def main():
  parser = argparse.ArgumentParser()
  parser.add_argument("infile")
  parser.add_argument("outdir")
  args = parser.parse_args()

  data, row_to_cell_line, gene_names = fl.parse_achilles(args.infile)
  for i, cell_line in enumerate(row_to_cell_line):
    data_cell_line = data[i,:]

    plt.clf()
    plt.hist(data_cell_line)
    plt.ylabel("Frequency")
    plt.xlabel("Probability")
    plt.title(cell_line)
    plt.savefig(os.path.join(args.outdir, "{}.png".format(cell_line)))

if __name__ == "__main__":
  main()
