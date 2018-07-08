#!/usr/bin/env python
import argparse
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import factorlib as fl

if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument('--raw')
  parser.add_argument('--U-out')
  parser.add_argument('--outfile')
  args = parser.parse_args()

  raw, row_names, col_names = fl.parse_achilles(args.raw)
  arr = np.genfromtxt(args.U_out, delimiter=",")
  rv = sns.clustermap(arr, yticklabels=row_names, z_score=0)
  plt.savefig(args.outfile, bbox_inches='tight')
