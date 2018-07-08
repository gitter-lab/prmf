#!/usr/bin/env python
# based on http://scikit-learn.org/stable/auto_examples/decomposition/plot_faces_decomposition.html#sphx-glr-auto-examples-decomposition-plot-faces-decomposition-py
# Authors: Aaron Baker, Vlad Niculae, Alexandre Gramfort
# License: BSD 3 clause
import os, os.path
import argparse
import numpy as np

import logging
from time import time

from numpy.random import RandomState
import matplotlib.pyplot as plt

from sklearn import decomposition
from sklearn.datasets import fetch_olivetti_faces
from sklearn.cluster import MiniBatchKMeans

import networkx as nx
import factorlib as fl

def main():
  parser = argparse.ArgumentParser(description="""
Run NMF following the scikit-learn demo on the Olivetti faces dataset.
Prepare AMPL files to solve the same problem.
""")
  parser.add_argument("outdir")
  args = parser.parse_args()

  comps_fp = os.path.join(args.outdir, "components.csv")
  data_image_fp = os.path.join(args.outdir, "olivetti_faces.png")
  comps_image_fp = os.path.join(args.outdir, "olivetti_comps.png")
  opt_fp = os.path.join(args.outdir, "optimize_log.txt")
  ampl_data_fp = os.path.join(args.outdir, "olivetti_data.dat")
  ampl_params_fp = os.path.join(args.outdir, "olivetti_params.dat")

  n_row, n_col = 2, 3
  n_components = n_row * n_col
  image_shape = (64, 64)
  rng = RandomState(0)
  # 400 (10 people x 40 images) x 4096 (64 x 64 images)
  dataset = fetch_olivetti_faces(shuffle=True, random_state=rng)
  faces = dataset.data
  n_samples, n_features = faces.shape

  # TODO transpose due to way my AMPL model is written (maybe should be rewritten to keep 
  # canonical obs x features data shape)
  with open(ampl_data_fp, 'w') as ampl_data_fh:
    fl.write_ampl_data(faces.transpose(), ampl_data_fh)
  with open(ampl_params_fp, 'w') as ampl_params_fh:
    fl.write_ampl_params(n_components, ampl_params_fh)

  # global centering
  faces_centered = faces - faces.mean(axis=0)
  # local centering
  faces_centered -= faces_centered.mean(axis=1).reshape(n_samples, -1)
  fl.plot_olivetti_components(plt, faces_centered[:n_components], "First centered Olivetti faces")
  plt.savefig(data_image_fp)

  # objective:
  # 0.5 * ||X - WH||_Fro^2
  #   + alpha * l1_ratio * ||vec(W)||_1
  #   + alpha * l1_ratio * ||vec(H)||_1
  #   + 0.5 * alpha * (1 - l1_ratio) * ||W||_Fro^2
  #   + 0.5 * alpha * (1 - l1_ratio) * ||H||_Fro^2
  name = 'Non-negative components - NMF'
  estimator = decomposition.NMF(n_components=n_components, init='nndsvda', tol=5e-3)
  center = False

  t0 = time()
  data = faces

  if center:
    data = faces_centered

  print("data shape {}".format(data.shape))
  estimator.fit_transform(data)
  train_time = (time() - t0)
  print("done in %0.3fs" % train_time)

  if hasattr(estimator, 'cluster_centers_'):
    components_ = estimator.cluster_centers_
  else:
    components_ = estimator.components_

  # Plot an image representing the pixelwise variance provided by the
  # estimator e.g its noise_variance_ attribute. The Eigenfaces estimator,
  # via the PCA decomposition, also provides a scalar noise_variance_
  # (the mean of pixelwise variance) that cannot be displayed as an image
  # so we skip it.
  #if (hasattr(estimator, 'noise_variance_') and
  #    estimator.noise_variance_.ndim > 0):  # Skip the Eigenfaces case
  #  plot_gallery("Pixelwise variance",
  #         estimator.noise_variance_.reshape(1, -1), n_col=1,
  #         n_row=1)
  fl.plot_olivetti_components(plt, components_[:n_components], '%s - Train time %.1fs' % (name, train_time))
  plt.savefig(comps_image_fp)

  # n_components x n_col
  np.savetxt(comps_fp, estimator.components_, delimiter=",")

  # also write components for AMPL
  # transform comps to n_obs x n_latent
  Gs = []
  comps = estimator.components_.transpose()
  n_obs, n_latent = comps.shape
  print("comps.shape:")
  print(comps.shape)
  for i in range(n_latent):
    latent_factor = comps[:,i]
    G = fl.vec_to_graph_prim(latent_factor)
    Gs.append(G)
  with open(ampl_data_fp, 'a') as ampl_data_fh:
    # TODO 
    # write_ampl_laplacians(Ls, ampl_data_fh)
    for i in range(len(Gs)):
      G = Gs[i]
      G_path = os.path.join(args.outdir, "manifold_{}.graphml".format(i))
      nx.write_graphml(G, G_path)

  # TODO write to log file with other optimization components
  opt_fh = open(opt_fp, 'w')
  opt_fh.write("reconstruction_err_: {}\n".format(estimator.reconstruction_err_))

if __name__ == "__main__":
  main()
