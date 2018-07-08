#!/usr/bin/env python
import argparse
import os, os.path
import math
import numpy as np
from scipy.stats import gamma
import matplotlib.pyplot as plt

def main():
  parser = argparse.ArgumentParser()
  parser.add_argument("outfile")
  args = parser.parse_args()

  ts = list(map(lambda x: x + 1, range(10)))
  plt.clf()
  xs = np.linspace(gamma.ppf(0.01, ts[0], 0, 1), gamma.ppf(0.99, ts[-1], 0, 1), 100)
  alpha_portion = 0.9
  for t in ts:
    # set location of pdf by t = alpha * beta
    alpha = math.pow(t, alpha_portion)
    beta = math.pow(t, 1-alpha_portion)
    ys = gamma.pdf(xs, alpha, 0, beta)
    plt.plot(xs, ys, label='t={}'.format(t))

  plt.title("Plot of X ~ Gamma(t^a, t^(1-a)) with a = 0.9")
  plt.xlabel("x")
  plt.ylabel("P(X = x)")
  plt.legend()
  plt.savefig(args.outfile, bbox_inches='tight')

if __name__ == "__main__":
  main()
