#!/usr/bin/env python
import sys, argparse
import numpy as np
import factorlib as fl
from sklearn.metrics import mean_squared_error as mse

def compute_mse(Z_true, Z_pred):
  nr, n_true = Z_true.shape
  nr, n_pred = Z_pred.shape
  mse_arr = np.zeros((n_true, n_pred))
  for i in range(n_true):
    for j in range(n_pred):
      mse_arr[i,j] = mse(Z_true[:,i], Z_pred[:,j])
  return mse_arr

def compute_corr(Z_true, Z_pred):
  nr, n_true = Z_true.shape
  nr, n_pred = Z_pred.shape
  corr_arr = np.zeros((n_true, n_pred))
  for i in range(n_true):
    for j in range(n_pred):
      corr_rv = np.corrcoef(Z_true[:,i], Z_pred[:,j])
      corr_arr[i,j] = corr_rv[0,1]
  return corr_arr

if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument('--true-z', required=True)
  parser.add_argument('--pred-z', required=True)
  parser.add_argument('--matching', '-m', action='store_true', help="""If set, perform a one-to-one 
matching between columns of --true-z and columns of --pred-z; otherwise, many columns of 
--pred-z may match against the same column of --true-z.""")
  parser.add_argument('--correlation', '-c', action='store_true', help="""If set, compute correlation
between pairs of column vectors. Otherwise, compute mean squared error.""")
  args = parser.parse_args()

  Z_true = np.genfromtxt(args.true_z, delimiter=",")
  Z_pred = np.genfromtxt(args.pred_z, delimiter=",")
  nr, n_pred = Z_pred.shape

  # true x pred
  perf_arr = None
  if args.correlation:
    perf_arr = compute_corr(Z_true, Z_pred)
  else:
    perf_arr = compute_mse(Z_true, Z_pred)

  # column indexes may be shuffled between Z_true and Z_pred
  perf_best = np.zeros((perf_arr.shape[1],))
  perf_best_ind = np.zeros((perf_arr.shape[1],))
  if args.matching:
    # do not allow columns of Z_true to be reused
    G, factor_node_ids, pathway_node_ids = fl.mat_to_bipartite(perf_arr)
    matching = nx.max_weight_matching(G)
    match_tpls = transform_matching(G, matching, factor_node_ids)
    for id1, id2, weight in match_tpls:
      true_ind = id1[1:]
      pred_ind = id2[1:]
      perf_best[pred_ind] = weight
      perf_best_ind[pred_ind] = true_ind
  else:
    # allow one column of Z_true to be reused
    perf_best = np.min(perf_arr, axis=0)
    perf_best_ind = np.argmin(perf_arr, axis=0)

  if not args.correlation:
    # normalize mse by column's 2 norm
    perf_best_normal = np.zeros(perf_best.shape)
    for i in range(n_pred):
      perf_best_normal[i] = perf_best[i] / np.linalg.norm(Z_true[:,perf_best_ind[i]])
    perf_best = perf_best_normal

  print('\t'.join(map(str, perf_best)))
  if args.correlation:
    # then summary metric is maximum
    print(np.max(perf_best))
  else:
    # then summary metric is mean
    print(np.mean(perf_best))
