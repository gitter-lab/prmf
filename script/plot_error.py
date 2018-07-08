#!/usr/bin/env python
import os, os.path
import sys, argparse
import matplotlib.pyplot as plt
import numpy as np

def parse_error(fh):
  rv = {}
  for line in fh:
    line = line.rstrip()
    words = line.split(",")
    ind0 = int(words[0])
    ind1 = int(words[1])
    error = float(words[2])
    rv[(ind0, ind1)] = error
  return rv

if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument("--error-init", "-i", required=True, type=argparse.FileType('r'))
  parser.add_argument("--error_m3r", "-m", required=True, type=argparse.FileType('r'))
  parser.add_argument("--outfile", "-o", required=True)
  args = parser.parse_args()

  init_map = parse_error(args.error_init)
  m3r_map = parse_error(args.error_m3r)

  keys = sorted(init_map.keys())
  xs = np.array(list(range(len(keys))))
  ys_init = np.array(list(map(lambda x: init_map[x], keys)))
  ys_m3r = np.array(list(map(lambda x: m3r_map[x], keys)))

  min_v = np.min(np.concatenate((ys_init, ys_m3r)))
  max_v = np.max(np.concatenate((ys_init, ys_m3r)))

  plt.xlim((min_v, max_v))
  plt.ylim((min_v, max_v))
  plt.scatter(ys_init, ys_m3r)
  plt.title('Pixel-specific error')
  plt.xlabel('nmf_init')
  plt.ylabel('nmf_m3r')
  plt.plot(xs, xs, 'k-')
  plt.savefig(args.outfile)
