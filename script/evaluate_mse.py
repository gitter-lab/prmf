#!/usr/bin/env python
import sys, argparse
import numpy as np
from sklearn.metrics import mean_squared_error as mse

if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument('--true-z', required=True)
  parser.add_argument('--pred-z', required=True)
  parser.add_argument('--outdir', required=True)
  args = parser.parse_args()

  Z_true = np.genfromtxt(args.true_z, delimiter=",")
  Z_pred = np.genfromtxt(args.pred_z, delimiter=",")

  # true x pred
  nr, n_true = Z_true.shape
  nr, n_pred = Z_pred.shape
  mse_arr = np.zeros((n_true, n_pred))
  for i in range(n_true):
    for j in range(n_pred):
      mse_arr[i,j] = mse(Z_true[:,i], Z_pred[:,j])

  # column indexes may be shuffled between Z_true and Z_pred
  # allow one column of Z_true to be reused
  mse_best = np.min(mse_arr, axis=0)
  mse_best_ind = np.argmin(mse_arr, axis=0)

  # normalize mse by column's 2 norm
  mse_best_normal = np.zeros(mse_best.shape)
  for i in range(n_pred):
    mse_best_normal[i] = mse_best[i] / np.linalg.norm(Z_true[mse_best_ind[i]])

  print('\t'.join(map(str, mse_best_normal)))
  print(np.mean(mse_best_normal))
