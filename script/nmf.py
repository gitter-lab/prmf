#!/usr/bin/env python
import sys, argparse
import os, os.path
import numpy as np
from sklearn import decomposition

def main():
  parser = argparse.ArgumentParser(description="""
Run traditional NMF:
  min_{U,V > 0} || X - UV^T ||_F

where X has shape (n_obs, n_feature)
  U has shape (n_obs, n_latent)
  V has shape (n_feature, n_latent)
""")
  parser.add_argument("--data", type=str, help="Input data matrix in CSV format", required=True)
  parser.add_argument("--outdir", "-o", type=str, help="Directory to write results to", required=True)
  parser.add_argument("--k-components", type=int, default=6)
  parser.add_argument("--seed", default=None, help="RNG seed")
  args = parser.parse_args()

  U_fp = os.path.join(args.outdir, "U.csv")
  V_fp = os.path.join(args.outdir, "V.csv")
  obj_fp = os.path.join(args.outdir, "obj.csv")

  X = np.genfromtxt(args.data, delimiter=",")

  # normalize assuming m < n (true for my application but unusual)
  m, n = X.shape
  if m > n:
    X = X.transpose()

  estimator = None
  if args.seed is None:
    estimator = decomposition.NMF(n_components=args.k_components)
  else:
    estimator = decomposition.NMF(n_components=args.k_components, random_state=int(args.seed))

  U = estimator.fit_transform(X)
  V = estimator.components_.astype(np.double).transpose()
  obj = estimator.reconstruction_err_
  #obj = np.linalg.norm(X - U * V_t)

  np.savetxt(U_fp, U, delimiter=",")
  np.savetxt(V_fp, V, delimiter=",")
  with open(obj_fp, 'w') as fh:
    fh.write("{} = {}".format("obj", obj))

if __name__ == "__main__":
  main()
