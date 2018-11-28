#!/usr/bin/env python
import argparse, sys
import factorlib as fl
import re
import os, os.path

def main():
  parser = argparse.ArgumentParser()
  parser.add_argument('--indir', help='Directory containing nmf_pathway.py outputs for multiple different runs, each with a different initialization', required=True)
  args = parser.parse_args()

  # <indir>/run\d+/obj.txt
  run_regexp = re.compile(r'run\d+')
  run_dirs = []
  for fname in os.listdir(args.indir):
    if run_regexp.match(fname) is not None:
      run_dirs.append(os.path.join(args.indir, fname))

  k = None
  counts = []
  for run_dir in run_dirs:
    converge = None
    init = None

    obj_path = os.path.join(run_dir, 'obj.txt')
    if os.path.exists(obj_path):
      converge = fl.parse_pathway_obj(obj_path)
    else:
      sys.stderr.write("[warning] Missing objective output for run with output directory {}\n".format(os.path.basename(obj_path)))

    init_path = os.path.join(run_dir, 'init_pathways.txt')
    if os.path.exists(init_path):
      init = fl.parse_init(init_path)
    else:
      sys.stderr.write("[warning] Missing initialization output for run with output directory {}\n".format(os.path.basename(obj_path)))

    if init is not None and converge is not None:
      # then include results in histogram
      # count number of init pathways which persist in the converged set
      count = 0
      converge_set = set(converge.items())
      for init_pathway in init:
        if init_pathway in converge_set:
          count += 1
      counts.append(count)

    if k is None:
      k = len(init)

  # make a histogram of counts
  n_bins = k+1
  n, bins, patches = plt.hist(count, bins=n_bins)

  # TODO
  plt.xlabel()
  plt.savefig()


if __name__ == "__main__":
  main()