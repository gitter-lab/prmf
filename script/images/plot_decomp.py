#!/usr/bin/env python
import argparse
import math
import copy
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from sklearn import decomposition

from matplotlib import rc
rc('text', usetex=True)

DPI = 96
LINE_WIDTH = 2.0
LABEL_FONT_SIZE = 'xx-large'

def plot_decomp(feature_vecs, aspect_ratio, ax_in, n_feature=5, n_latent=4):
  """
  Parameters
  ----------
  feature_vecs : list of (np.array or None)
    a list of size <n_latent> which specifies the latent vectors of the matrix W
  """
  feature_vecs = copy.copy(feature_vecs)

  n_obs = math.ceil(n_feature * aspect_ratio)

  total_width = n_feature + n_latent + n_feature
  approx_portion = total_width / 12
  total_width += approx_portion
  feature_fraction = n_feature / total_width
  approx_fraction = approx_portion / total_width
  latent_fraction =  n_latent / total_width
  latent_fraction_h = latent_fraction * aspect_ratio

  left, bottom, right, top = ax_in.get_position().get_points().flatten()
  gs = mpl.gridspec.GridSpec(2, 4, 
    left=left, bottom=bottom, right=right, top=top,
    width_ratios=[feature_fraction, approx_fraction, latent_fraction, feature_fraction],
    height_ratios=[latent_fraction_h, 1 - latent_fraction_h])

  ax = plt.subplot(gs[:,0])
  data = np.random.rand(n_obs, n_feature)
  for i, feature_vec in enumerate(feature_vecs):
    if feature_vec is None:
      feature_vecs[i] = np.random.rand(n_feature, 1)

  # X = WH
  # X is shape n_obs, n_features
  # W is shape n_obs, n_latent
  # H is shape n_latent, n_features
  H = np.concatenate(feature_vecs, axis=1).transpose()
  W = np.random.rand(n_obs, n_latent)
  X = np.dot(W, H)

  mappable = ax.pcolor(X, cmap=plt.cm.Reds, linewidth=LINE_WIDTH)
  ax.set_title("Features", fontsize=LABEL_FONT_SIZE)
  ax.set_ylabel("Observations", fontsize=LABEL_FONT_SIZE)
  ax.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
  ax.tick_params(axis='y', which='both', left=False, right=False, labelleft=False)

  ax = plt.subplot(gs[:,1])
  ax.axis('off')
  ax.text(0.0, 0.5, r"$\approx$", fontsize=60)

  ax = plt.subplot(gs[:,2])
  mappable = ax.pcolor(W, cmap=plt.cm.Reds, linewidth=LINE_WIDTH)
  ax.set_title("Latent", fontsize=LABEL_FONT_SIZE)
  #ax.set_ylabel("Observations", fontsize=LABEL_FONT_SIZE)
  ax.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
  ax.tick_params(axis='y', which='both', left=False, right=False, labelleft=False)

  ax = plt.subplot(gs[0,3])
  mappable = ax.pcolor(H, cmap=plt.cm.Reds, linewidth=LINE_WIDTH)
  ax.set_title("Features", fontsize=LABEL_FONT_SIZE)
  #ax.set_ylabel("Latent", fontsize=LABEL_FONT_SIZE)
  ax.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
  ax.tick_params(axis='y', which='both', left=False, right=False, labelleft=False)

def main():
  parser = argparse.ArgumentParser()
  parser.add_argument("out")
  args = parser.parse_args()

  fig_width = 900
  fig_height = 600
  aspect_ratio = fig_width / fig_height
  fig = plt.figure(figsize=(fig_width/DPI, fig_height/DPI), dpi=DPI)
  ax = fig.gca()

  n_feature = 5 
  n_latent = 4 
  feature_vecs = []
  for i in range(n_latent):
    data = np.random.rand(n_feature, 1)
    feature_vecs.append(data)
  plot_decomp(feature_vecs, aspect_ratio, ax, n_feature=n_feature, n_latent=n_latent)

  fig.savefig(args.out)

if __name__ == "__main__":
  main()
