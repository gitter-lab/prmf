#!/usr/bin/env python
import argparse, sys
import os, os.path
import re
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

def parse_eval(fp):
  rv = None
  with open(fp, 'r') as fh:
    i = 0
    for line in fh:
      line = line.rstrip()
      i += 1
      if i == 2:
        rv = float(line)
  return rv

if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument('--indir', '-i', help='Directory of results to plot')
  parser.add_argument('--outdir', '-o', help='Directory to place plots')
  args = parser.parse_args()

  # run/method/evaluate_mse.out
  run_regexp = re.compile(r'run\d+')
  nmf_rel_path = os.path.join("NMF", "evaluate_mse.out")
  plier_rel_path = os.path.join("PLIER", "evaluate_mse.out")

  xs = []
  ys = []
  for fname in os.listdir(args.indir):
    match_data = run_regexp.match(fname)
    if match_data is not None:
      nmf_eval_fp = os.path.join(args.indir, fname, nmf_rel_path)
      plier_eval_fp = os.path.join(args.indir, fname, plier_rel_path)

      nmf_eval_v = parse_eval(nmf_eval_fp)
      plier_eval_v = parse_eval(plier_eval_fp)

      if not nmf_eval_v is None and not plier_eval_v is None:
        xs.append(nmf_eval_v)
        ys.append(plier_eval_v)
  xs = np.array(xs)
  ys = np.array(ys)

  x_max = np.percentile(xs, 95)
  y_max = np.percentile(ys, 95)
  both_max = np.max([x_max, y_max])

  plt.scatter(xs, ys, linewidths=2.0)
  plt.plot(np.linspace(0, both_max), np.linspace(0, both_max), 'k-')
  plt.xlim([0, both_max])
  plt.ylim([0, both_max])
  plt.xlabel('NMF')
  plt.ylabel('PLIER')
  plt.title('Average MSE from PLIER-based Simulation')

  ofp = os.path.join(args.outdir, 'mse_plot.png')
  plt.savefig(ofp)
