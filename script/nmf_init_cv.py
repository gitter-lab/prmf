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
import prmf

def parse_h(fp, delimiter="\t"):
  """
  Returns
  -------
  H : np.array
    array with smaller dimension along the row axis
  """
  H = np.genfromtxt(fp, delimiter=delimiter)
  n_row, n_col = H.shape
  if(n_col < n_row):
    H = H.transpose()
  return H

def perturb_h(H):
  # uniform random in [0,1)
  noise = np.random.rand(*H.shape)
  return H + noise

def main():
  parser = argparse.ArgumentParser(description="""
Write optimization performance to stdout and optimal W and H matrices to separate CSV files
""")
  parser.add_argument("--prior", type=str, default="prim", help="One of prim, comb")
  parser.add_argument("--k-components", type=int, default=6)
  parser.add_argument("--outdir", "-o", type=str, help="Directory to write results to")
  args = parser.parse_args()

  opt_out_fh = open(os.path.join(args.outdir, "opt.txt"), 'w')
  out_w_fh = open(os.path.join(args.outdir, "W_opt_init.csv"), 'wb')
  out_h_fh = open(os.path.join(args.outdir, "H_opt_init.csv"), 'wb')
  out_indexes_fh = open(os.path.join(args.outdir, 'indexes.csv'), 'w')
  prior_fh = open(os.path.join(args.outdir, 'H_prior.csv'), 'wb')
  true_fh = open(os.path.join(args.outdir, 'H_true.csv'), 'wb')
  pixel_error_fh = open(os.path.join(args.outdir, "error_init.csv"), 'w')
  test_fh = open(os.path.join(args.outdir, "X_test.csv"), "wb")

  # 400 (10 people x 40 images) x 4096 (64 x 64 images)
  rng = RandomState(0)
  dataset = fetch_olivetti_faces(shuffle=True, random_state=rng)
  data = dataset.data
  n_samples, n_features = data.shape

  # 50-50 cross validation
  X_train, X_test = sklearn.model_selection.train_test_split(data, train_size=0.5, stratify=dataset.target)
  np.savetxt(test_fh, X_test, delimiter=",")
  estimator = decomposition.NMF(n_components=args.k_components)
  W_init = estimator.fit(X_train) # this will actually be discarded
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
  prmf.write_ampl_data(X_test.transpose(), ampl_data_fh)
  prmf.write_ampl_params(n_components, ampl_params_fh)

  Gs = []
  for k in range(args.k_components):
    comp = H_init[k,:]
    G = None
    if args.prior == "prim":
      G = prmf.vec_to_graph_prim(comp)
    elif args.prior == "comb":
      G = prmf.vec_to_graph_comb(comp)
    else:
      sys.stderr.write("invalid --prior, using prim\n")
      G = prmf.vec_to_graph_prim(comp)
    Gs.append(G)
  laplacians = map(nx.laplacian_matrix, Gs)
  prmf.ampl_write_sparse_arrs(laplacians, ampl_lapl_fh)

  # }}--

if __name__ == "__main__":
  main()
