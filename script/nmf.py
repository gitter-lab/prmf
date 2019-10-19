#!/usr/bin/env python
import sys, argparse
import os, os.path
import numpy as np
import factorlib as fl
from sklearn import decomposition
import csv
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.model_selection import KFold

def main():
  parser = argparse.ArgumentParser(description="""
Run traditional NMF:
  min_{U,V > 0} || X - UV^T ||_F

where X has shape (n_obs, n_feature)
  U has shape (n_obs, n_latent)
  V has shape (n_feature, n_latent)

Writes U.csv with the U matrix including row and column names, V.csv with the V matrix including row and column names, and obj.csv which includes information about the objective function and cross validation if --cross-validation argument is provided.
""")
  parser.add_argument("--data", type=str, help="Input data matrix in CSV format", required=True)
  parser.add_argument("--outdir", "-o", type=str, help="Directory to write results to", required=True)
  parser.add_argument("--k-latent", type=int, default=6)
  parser.add_argument("--seed", default=None, help="RNG seed")
  parser.add_argument("--delimiter", "-d", default=",", type=str, help="Field delimiter used in the --data file")
  parser.add_argument("--cross-validation", "-c", type=float, help="If provided, use the --cross-validation value as a fraction of the samples to hold out and measure model performance with")
  args = parser.parse_args()

  U_fp = os.path.join(args.outdir, "U.csv")
  V_fp = os.path.join(args.outdir, "V.csv")
  obj_fp = os.path.join(args.outdir, "obj.csv")

  # assume csv has row and column names
  X = pd.read_csv(args.data, sep=args.delimiter, header='infer', index_col=0)
  samples = list(X.index)
  nodelist = list(X.columns)
  X = X.values

  # normalize assuming m < n (true for my application but unusual)
  m, n = X.shape
  if m > n:
    X = X.transpose()

  estimator = None
  if args.seed is None:
    estimator = decomposition.NMF(n_components=args.k_latent, init='random')
  else:
    estimator = decomposition.NMF(n_components=args.k_latent, random_state=int(args.seed), init='random')

  # cross validation
  # TODO validate cross_validation input in [0,1]
  # TODO change cross_validation input for KFold
  X_test = None
  X_train = None
  if args.cross_validation is not None:
    kf = KFold(n_splits = round(1/args.cross_validation))
    for train_index, test_index in kf.split(X):
      X_train = X[train_index]
      X_test = X[test_index]

      X = X_train
      samples = [samples[i] for i in train_index]
      break

  U = estimator.fit_transform(X)
  V = estimator.components_.astype(np.double).transpose()
  obj = estimator.reconstruction_err_
  #obj = np.linalg.norm(X - U * V_t)

  # cross validation
  avg_normalized_test_error = None
  if args.cross_validation is not None:
    normalized_test_errors = fl.measure_cv_performance(V, X_test)
    avg_normalized_test_error = np.mean(normalized_test_errors)
    error_fp = os.path.join(args.outdir, 'test_error.csv')
    np.savetxt(error_fp, normalized_test_errors, delimiter=",")

  np.savetxt(U_fp, U, delimiter=",")
  np.savetxt(V_fp, V, delimiter=",")

  U = pd.DataFrame(U, index=samples, columns=list(map(lambda x: "LV{}".format(x), range(args.k_latent))))
  V = pd.DataFrame(V, index=nodelist, columns=list(map(lambda x: "LV{}".format(x), range(args.k_latent))))
  U.to_csv(U_fp, sep=",", index=True, quoting=csv.QUOTE_NONNUMERIC)
  V.to_csv(V_fp, sep=",", quoting=csv.QUOTE_NONNUMERIC)

  with open(obj_fp, 'w') as fh:
    fh.write("{} = {}\n".format("obj", obj))
    if avg_normalized_test_error is not None:
      fh.write("{} = {}".format("avg_normalized_test_error", avg_normalized_test_error))

if __name__ == "__main__":
  main()
