#!/usr/bin/env python
import sys, argparse
import os, os.path
import numpy as np
import scipy.stats as st
import random
from sklearn import decomposition

def main():
  parser = argparse.ArgumentParser(description="""
Run traditional NMF:
  min_{U,V > 0} || X - UV^T ||_F

where X has shape (n_obs, n_feature)
  U has shape (n_obs, k_components)
  V has shape (n_feature, k_components)

on simulated data X = YZ^T where
  Y ~ Gamma(5,1)
  Z = A / n
  A ~ Multinomial(n_feature, (1/n_feature, 1/n_feature, ...))
""")
  parser.add_argument("--outdir", "-o", type=str, help="Directory to write results to", required=True)
  parser.add_argument("--n-obs", type=int, default=300)
  parser.add_argument("--k-components", type=int, default=6)
  parser.add_argument("--n-features", type=int, default=20000)
  parser.add_argument("--precision", type=int, default=8, help="Number of decimal points in proportions in Z")
  parser.add_argument("--gamma-shape", type=float, default=5, help="Gamma shape parameter")
  parser.add_argument("--init", type=bool, default=False, help="Use custom initialization")
  parser.add_argument("--seed", default=None, help="RNG seed")
  args = parser.parse_args()

  Y = np.zeros((args.n_obs, args.k_components))
  for k in range(args.k_components):
    Y[:,k] = st.gamma.rvs(args.gamma_shape, size=args.n_obs)
  Z = st.multinomial.rvs(args.precision, np.repeat(1/args.n_features, args.n_features), size=args.k_components)
  # add pseudocounts to Z
  Z = Z + 1
  total = args.precision + args.n_features
  Z = Z / total
  X_true = np.dot(Y,Z)

  estimator = None
  init = 'nndsvdar'
  if(args.init):
    init = 'custom'
  if args.seed is None:
    estimator = decomposition.NMF(n_components=args.k_components, init=init)
  else:
    estimator = decomposition.NMF(n_components=args.k_components, random_state=int(args.seed), init=init)

  U = None
  if(args.init):
    # test prior-based initialization using simulated data
    # W_init follows Wu strategy
    # H_init gets ground truth
    # NOTE so sklearn.decomposition has X - WH
    # W := U
    # H := V^T
    population = list(range(args.n_features))
    sample = random.sample(population, args.k_components)
    W_init = X_true[:, sample].copy(order='C').astype(np.double)
    H_init = Z

    U = estimator.fit_transform(X_true, W=W_init, H=H_init)
  else:
    U = estimator.fit_transform(X_true)

  V = estimator.components_.astype(np.double).transpose()
  obj = estimator.reconstruction_err_
  # NOTE obj := np.linalg.norm(X_true - np.dot(U, V.transpose()))
  print(obj / np.linalg.norm(X_true))

if __name__ == "__main__":
  main()
