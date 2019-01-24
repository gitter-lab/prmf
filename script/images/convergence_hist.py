#!/usr/bin/env python
import argparse, sys
import re
import os, os.path
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import factorlib as fl

def main():
  parser = argparse.ArgumentParser(description="""
Create a histogram of converged pathways
""")
  parser.add_argument('--indir', help='Directory containing nmf_pathway.py outputs for multiple different runs, each with a different initialization', required=True)
  parser.add_argument('--outdir', help='Directory to place image file', required=True)
  args = parser.parse_args()

  # <indir>/run\d+/obj.txt
  run_regexp = re.compile(r'run\d+')
  run_dirs = []
  for fname in os.listdir(args.indir):
    if run_regexp.match(fname) is not None:
      run_dirs.append(os.path.join(args.indir, fname))

  N = len(run_dirs)

  # count number of converged pathways
  pathway_to_count = {}
  for run_dir in run_dirs:
    pathways_dict = fl.parse_pathway_obj(os.path.join(run_dir, 'obj.txt'))
    pathways = pathways_dict.values()
    pathway_basenames = list(map(os.path.basename, pathways))
    for pathway in pathway_basenames:
      if pathway in pathway_to_count:
        pathway_to_count[pathway] += 1 
      else:
        pathway_to_count[pathway] = 1
  pathway_count_pairs = sorted(pathway_to_count.items(), key=lambda x: x[1], reverse=True)
  hist_csv_fp = os.path.join(args.outdir, 'convergence_hist.csv')
  with open(hist_csv_fp, 'w') as hist_csv_fh:
    for pathway, count in pathway_count_pairs:
      hist_csv_fh.write("{},{}\n".format(pathway, count))

  pathway_to_freq = {}
  for pathway in pathway_to_count.keys():
    pathway_to_freq[pathway] = pathway_to_count[pathway] / N

  # histogram for converged pathways - {{
  freqs = pathway_to_freq.values()
  freqs_per = np.array(list(freqs)) * 100

  plt.clf()
  n, bins, patches = plt.hist(freqs_per, bins=np.arange(0,110,10))
  hist_fp = os.path.join(args.outdir, 'convergence_hist.png')
  plt.xlabel('Pathway Frequency')
  plt.ylabel('Frequency')
  plt.title('Stability of Converged Pathways')
  plt.savefig(hist_fp)
  # }} -

if __name__ == "__main__":
  main()
