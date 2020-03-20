#!/usr/bin/env python
import sys, argparse
import os, os.path
import re
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

RUN_REGEXP = re.compile(r'run\d+')
OBJ_AVP_REGEXP = re.compile('(\w+)\W+(\w+)')
ATTRIBUTES = ['obj', 'fro', 'manifold', 'delta', 'ignore', 'gamma', 'recon']

def get_stats(run_dir):
  obj_fp = os.path.join(run_dir, 'obj.txt')

  # first 7 lines contain attribute = value pairs
  obj_func_parts = {}
  with open(obj_fp) as fh:
    line_no = 0
    for line in fh:
      line = line.rstrip()
      line_no += 1
      if line_no <= 7:
        match_data = OBJ_AVP_REGEXP.search(line)
        if match_data is not None:
          attribute = match_data.group(1)
          value = match_data.group(2)
          obj_func_parts[attribute] = float(value)

  return obj_func_parts

def main():
  parser = argparse.ArgumentParser(description="""
Reduce multiple runs of nmf_pathway.py (via nmf_pathway_rr.py)

TODO
----
Currently this script just collects statistics about the runs
""")
  parser.add_argument('--indir')
  parser.add_argument('--outdir')
  parser.add_argument('--out-format', default='hist', help='If "hist" create a histogram for each objective function attribute')
  args = parser.parse_args()

  run_dirs = []
  for fname in os.listdir(args.indir):
    fpath = os.path.join(args.indir, fname)
    match_data = RUN_REGEXP.search(fpath)
    if match_data is not None:
      run_dirs.append(fpath)

  all_stats = []
  for run_dir in run_dirs:
    stats = get_stats(run_dir)
    all_stats.append(stats)

  if args.out_format == 'hist':
    for attribute in ATTRIBUTES:
      ofp = os.path.join(args.outdir, '{}.png'.format(attribute))

      counts = []
      for stats in all_stats:
        if attribute in stats:
          counts.append(stats[attribute])
        else:
          sys.stderr.write('Invalid attribute {} is not present in objective file output for run {}\n'.format(attribute, run_no))
          sys.exit(22)
      counts = np.array(counts)

      n_bins = 20
      plt.clf()
      plt.title('Frequency of {}'.format(attribute))
      n_row = 3
      n_col = 3
      fig, axes = plt.subplots(n_row, n_col)
      axes = axes.flatten()
      master_ax = fig.gca()

      percentile_edges = np.linspace(0,100,num=10)
      for i in range(n_row * n_col):
        ax = axes[i]
        ax.hist(counts, range=(np.percentile(counts, percentile_edges[i]), np.percentile(counts, percentile_edges[i+1])))

      fig.savefig(ofp)

  else:
    sys.stderr.write('Invalid --out-format = {}\n'.format(args.out_format))
    sys.exit(21)

if __name__ == "__main__":
  main()
