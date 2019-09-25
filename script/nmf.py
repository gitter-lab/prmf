#!/usr/bin/env python
import sys, argparse
import os, os.path
import numpy as np
import factorlib as fl
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
  parser.add_argument("--k-latent", type=int, default=6)
  parser.add_argument("--seed", default=None, help="RNG seed")
  parser.add_argument("--cross-validation", "-c", type=float, help="If provided, use the --cross-validation value as a fraction of the samples to hold out and measure model performance with")
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
    estimator = decomposition.NMF(n_components=args.k_latent)
  else:
    estimator = decomposition.NMF(n_components=args.k_latent, random_state=int(args.seed))

  # cross validation
  X_test = None
  if args.cross_validation is not None:
    # TODO validate cross_validation input in [0,1]
    X_train, X_test = train_test_split(X, test_size=args.cross_validation)
    X = X_train

  U = estimator.fit_transform(X)
  V = estimator.components_.astype(np.double).transpose()
  obj = estimator.reconstruction_err_
  #obj = np.linalg.norm(X - U * V_t)

  # cross validation
  if args.cross_validation is not None:
    normalized_test_errors = fl.measure_cv_performance(V, X_test)
    avg_normalized_test_error = np.mean(normalized_test_errors)
    error_fp = os.path.join(args.outdir, 'test_error.csv')
    np.savetxt(error_fp, normalized_test_errors, delimiter=",")
    obj_data['average_normalized_test_error'] = avg_normalized_test_error

  np.savetxt(U_fp, U, delimiter=",")
  np.savetxt(V_fp, V, delimiter=",")
  with open(obj_fp, 'w') as fh:
    fh.write("{} = {}".format("obj", obj))

if __name__ == "__main__":
  main()
