#!/usr/bin/env python
import os, os.path
import sys, argparse
import numpy as np

if __name__ == "__main__":
  parser = argparse.ArgumentParser(description="""
Evaluate performance of m3r method. Idea is to continue where nmf_init.py left off.
""")
  parser.add_argument("indir")
  parser.add_argument("outdir")
  args = parser.parse_args()

  # TODO names are confusing; ampl calls n_feature x n_comp matrix W
  pred_fh = open(os.path.join(args.indir, "W_opt_m3r.csv"), 'rb')
  true_fh = open(os.path.join(args.indir, 'H_true.csv'), 'rb')
  indexes_fh = open(os.path.join(args.indir, 'indexes.csv'), 'r')
  pixel_error_fh = open(os.path.join(args.outdir, 'error_m3r.csv'), 'w')
  opt_out_fh = open(os.path.join(args.outdir, 'opt_m3r.txt'), 'w')

  # normalize H's to comp x feature
  H_true = np.genfromtxt(true_fh, delimiter=",")
  H = np.genfromtxt(pred_fh, delimiter=',').transpose()
  H_valid = []
  for line in indexes_fh:
    line = line.rstrip()
    H_valid.append(tuple(map(int, line.split(','))))

  # evaluate performance
  mse = 0.0
  n_elem = len(H_valid)
  for (i,j) in H_valid:
    true_v = H_true[i,j]
    pred_v = H[i,j]
    error = (true_v - pred_v) ** 2 
    pixel_error_fh.write(",".join(map(str, (i,j,error))) + "\n")
    mse += error
  mse = mse / n_elem
  opt_out_fh.write("mse = {}\n".format(mse))
