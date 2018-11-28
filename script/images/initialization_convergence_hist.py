#!/usr/bin/env python
import argparse, sys
import re
import os, os.path
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import factorlib as fl

# TODO basename is used for ground truth in case file path changes (in case where programs 
# are run on different machines)

def count_initialization_and_convergence_overlap(run_dirs):
  """
  For each run in <run_dirs>, count how many of the pathways used for initialization
  exist in the set of pathways that the method converged to
  """
  k = None
  counts = []
  for run_dir in run_dirs:
    converge = None
    init = None

    obj_path = os.path.join(run_dir, 'obj.txt')
    if os.path.exists(obj_path):
      converge = fl.parse_pathway_obj(obj_path)
    else:
      sys.stderr.write("[warning] Missing objective output for run with output directory {}\n".format(os.path.basename(run_dir)))

    init_path = os.path.join(run_dir, 'init_pathways.txt')
    if os.path.exists(init_path):
      init = fl.parse_init(init_path)
    else:
      sys.stderr.write("[warning] Missing initialization output for run with output directory {}\n".format(os.path.basename(run_dir)))

    if init is not None and converge is not None:
      # then include results in histogram
      # count number of init pathways which persist in the converged set
      count = 0
      converge_set = set(converge.values())
      for init_pathway in init:
        if init_pathway in converge_set:
          count += 1
      counts.append(count)

      if k is None:
        k = len(init)
  return counts, k

def count_initialization_and_truth_overlap(run_dirs, ground_truth_fp):
  counts = []

  truth_set = set()
  with open(ground_truth_fp, 'r') as fh:
    for line in fh:
      line = line.rstrip()
      truth_set.add(os.path.basename(line))

  for run_dir in run_dirs:
    init = None

    init_path = os.path.join(run_dir, 'init_pathways.txt')
    if os.path.exists(init_path):
      init = fl.parse_init(init_path)
    else:
      sys.stderr.write("[warning] Missing initialization output for run with output directory {}\n".format(os.path.basename(run_dir)))

    if init is not None:
      # then include results in histogram
      # count number of init pathways which exist in the ground truth
      count = 0
      for init_pathway in init:
        if os.path.basename(init_pathway) in truth_set:
          count += 1
      counts.append(count)
  return counts

def count_converge_and_truth_overlap(run_dirs, ground_truth_fp):
  counts = []

  truth_set = set()
  with open(ground_truth_fp, 'r') as fh:
    for line in fh:
      line = line.rstrip()
      truth_set.add(os.path.basename(line))

  for run_dir in run_dirs:
    converge = None

    obj_path = os.path.join(run_dir, 'obj.txt')
    if os.path.exists(obj_path):
      converge = fl.parse_pathway_obj(obj_path)
    else:
      sys.stderr.write("[warning] Missing objective output for run with output directory {}\n".format(os.path.basename(run_dir)))

    if converge is not None:
      # then include results in histogram
      # count number of converged pathways which exist in the ground truth
      count = 0
      for converge_pathway in converge.values():
        if os.path.basename(converge_pathway) in truth_set:
          count += 1
      counts.append(count)
  return counts

def main():
  parser = argparse.ArgumentParser()
  parser.add_argument('--indir', help='Directory containing nmf_pathway.py outputs for multiple different runs, each with a different initialization', required=True)
  parser.add_argument('--ground-truth', help='File containing filepaths to the pathway files used as the basis of a simulation')
  parser.add_argument('--outdir', help='Directory to place image files', required=True)
  args = parser.parse_args()

  # <indir>/run\d+/obj.txt
  run_regexp = re.compile(r'run\d+')
  run_dirs = []
  for fname in os.listdir(args.indir):
    if run_regexp.match(fname) is not None:
      run_dirs.append(os.path.join(args.indir, fname))

  N = len(run_dirs)

  # histogram for initialization/convergence - {{
  # make a histogram of counts
  counts, k = count_initialization_and_convergence_overlap(run_dirs)

  plt.clf()
  n_bins = k+1
  n, bins, patches = plt.hist(counts, bins=np.arange(n_bins+1)-0.5)
  initial_converged_fp = os.path.join(args.outdir, 'initialization_convergence_hist.png')
  plt.xlabel('Number persisting')
  plt.ylabel('Frequency')
  plt.title('Initial pathways persisting in converged set (N = {})'.format(N))
  plt.savefig(initial_converged_fp)
  # }} -

  # histogram for initialization/ground truth - {{
  if args.ground_truth is not None:
    counts = count_initialization_and_truth_overlap(run_dirs, args.ground_truth)

    plt.clf()
    n_bins = k+1
    n, bins, patches = plt.hist(counts, bins=np.arange(n_bins+1)-0.5)
    initial_converged_fp = os.path.join(args.outdir, 'initialization_truth_hist.png')
    plt.xlabel('Number existing')
    plt.ylabel('Frequency')
    plt.title('Initial pathways existing in ground truth set (N = {})'.format(N))
    plt.savefig(initial_converged_fp)
  # }} -

  # histogram for ground truth/convergence - {{
  if args.ground_truth is not None:
    counts = count_converge_and_truth_overlap(run_dirs, args.ground_truth)

    plt.clf()
    n_bins = k+1
    n, bins, patches = plt.hist(counts, bins=np.arange(n_bins+1)-0.5)
    initial_converged_fp = os.path.join(args.outdir, 'convergence_truth_hist.png')
    plt.xlabel('Number existing')
    plt.ylabel('Frequency')
    plt.title('Converged pathways existing in ground truth set (N = {})'.format(N))
    plt.savefig(initial_converged_fp)
  # }} -

if __name__ == "__main__":
  main()
