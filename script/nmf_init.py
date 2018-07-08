#!/usr/bin/env python
import os, os.path
import sys, argparse
import random
import itertools as it
import numpy as np
from numpy.random import RandomState
from sklearn import decomposition
from sklearn.datasets import fetch_olivetti_faces
import sklearn.model_selection
import networkx as nx
import factorlib as fl

def perturb_h(H):
  # uniform random in [0,1)
  noise = np.random.rand(*H.shape)
  return H + noise

def main():
  parser = argparse.ArgumentParser(description="""
Run NMF with prior-based initialization
""")
  parser.add_argument("--data", type=str, help="Input data matrix in CSV format", required=True)
  parser.add_argument("--prior", type=str, help="feature-feature graph", required=True)
  parser.add_argument("--outdir", "-o", type=str, help="Directory to write results to", required=True)
  parser.add_argument("--k-components", type=int, default=6)
  parser.add_argument("--seed", default=None, help="RNG seed")
  args = parser.parse_args()

  if args.seed is not None:
    rng = RandomState(int(args.seed))

  opt_out_fh = open(os.path.join(args.outdir, "opt.txt"), 'w')
  out_w_fh = open(os.path.join(args.outdir, "W_opt_init.csv"), 'wb')
  out_h_fh = open(os.path.join(args.outdir, "H_opt_init.csv"), 'wb')
  out_indexes_fh = open(os.path.join(args.outdir, 'indexes.csv'), 'w')
  prior_fh = open(os.path.join(args.outdir, 'H_prior.csv'), 'wb')
  true_fh = open(os.path.join(args.outdir, 'H_true.csv'), 'wb')
  pixel_error_fh = open(os.path.join(args.outdir, "error_init.csv"), 'w')
  test_fh = open(os.path.join(args.outdir, "X_test.csv"), "wb")

  data = np.genfromtxt(args.data, delimiter=",")

  estimator = decomposition.NMF(n_components=args.k_components)
  W_init = estimator.fit(data)
  H_true = estimator.components_.astype(np.double)
  H_init = perturb_h(H_true)
  np.savetxt(prior_fh, H_init, delimiter=",")
  np.savetxt(true_fh, H_true, delimiter=",")

  # select 1 / 5 of prior, "validation set", to zero out and use to evaluate performance
  n_components, n_features = H_init.shape
  population = list(it.product(range(n_components), range(n_features)))
  H_valid = random.sample(population, int(0.2 * n_components * n_features))
  for (i,j) in H_valid:
    H_init[i,j] = 0
    out_indexes_fh.write(",".join(map(str,(i,j))) + "\n")

  # --{{ begin baseline 

  # sample n_components rows from data to form H
  # ||X - WH||
  # X = n_samples, n_features
  # W = n_samples, n_components
  # H = n_components, n_features
  population = list(range(n_samples))
  sample = random.sample(population, n_components)
  W_init = X_test[:, sample].copy(order='C').astype(np.double)

  estimator = decomposition.NMF(n_components=n_components, init='custom')
  W = estimator.fit_transform(X_test, W=W_init, H=H_init)
  H = estimator.components_

  # write other aspects of the model
  residual = estimator.reconstruction_err_
  opt_out_fh.write("residual = {}\n".format(residual))
  np.savetxt(out_w_fh, W, delimiter=",")
  np.savetxt(out_h_fh, H, delimiter=",")

  # evaluate performance
  mse = 0.0
  n_elem = len(H_valid)
  for (i,j) in H_valid:
    true_v = H_true[i,j]
    pred_v = H[i,j]
    error = (true_v - pred_v) ** 2 
    pixel_error_fh.write(",".join(map(str, (i,j,error))) + "\n")
    mse += error
  mse = mse / n_elem
  opt_out_fh.write("mse = {}\n".format(mse))

  # }}-- end baseline

  # --{{ begin multi manifold min regularization method

  # prepare files for ampl: olivetti.dat, laplacians.dat, params.dat
  ampl_data_fh = open(os.path.join(args.outdir, "olivetti.dat"), 'w')
  ampl_params_fh = open(os.path.join(args.outdir, "params.dat"), 'w')
  ampl_lapl_fh = open(os.path.join(args.outdir, "laplacians.dat"), 'w')

  # olivetti.dat writes X_test as n_feature x n_samples
  fl.write_ampl_data(X_test.transpose(), ampl_data_fh)
  fl.write_ampl_params(n_components, ampl_params_fh)

  Gs = []
  for k in range(args.k_components):
    comp = H_init[k,:]
    G = None
    if args.prior == "prim":
      G = fl.vec_to_graph_prim(comp)
    elif args.prior == "comb":
      G = fl.vec_to_graph_comb(comp)
    else:
      sys.stderr.write("invalid --prior, using prim\n")
      G = fl.vec_to_graph_prim(comp)
    Gs.append(G)
  laplacians = map(nx.laplacian_matrix, Gs)
  fl.ampl_write_sparse_arrs(laplacians, ampl_lapl_fh)

  # }}--

if __name__ == "__main__":
  main()
