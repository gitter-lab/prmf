#!/usr/bin/env python
import argparse
import copy
import math
import random
import networkx as nx
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

# https://www.infobyip.com/detectmonitordpi.php
DPI = 96
LINE_WIDTH = 2.0
LABEL_FONT_SIZE = 'xx-large'
TICK_FONT_SIZE = 'x-large'
mpl.rcParams['lines.markeredgewidth'] = 5.0
#'scatter.marker'

# TODO this is a copy from plot_decomp
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
  ax.text(-0.4, 0.5, r"$\approx$", fontsize=60)

  ax = plt.subplot(gs[:,2])
  mappable = ax.pcolor(W, cmap=plt.cm.Reds, linewidth=LINE_WIDTH)
  ax.set_title("Latent", fontsize=LABEL_FONT_SIZE)
  #ax.set_ylabel("Observations", fontsize=LABEL_FONT_SIZE)
  ax.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
  ax.tick_params(axis='y', which='both', left=False, right=False, labelleft=False)

  ax = plt.subplot(gs[0,3])
  ax.invert_yaxis()
  mappable = ax.pcolor(H, cmap=plt.cm.Reds, linewidth=LINE_WIDTH)
  ax.set_title("Features", fontsize=LABEL_FONT_SIZE)
  #ax.set_ylabel("Latent", fontsize=LABEL_FONT_SIZE)
  ax.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
  ax.tick_params(axis='y', which='both', left=False, right=False, labelleft=False)

def plot_latent_and_graph(G, fig, ax, column_index, padding=0.02, left_percent=0.3, pos=None, axis='off'):
  right_percent = 1 - left_percent

  ax.get_xaxis().set_visible(False)
  ax.get_yaxis().set_visible(False)

  # create labels centered at each heatmap block
  x0, y0, x1, y1 = ax.get_position().get_points().flatten()
  width = x1 - x0
  height = y1 - y0
  ax_1 = fig.add_axes([x0, y0, width * left_percent, height])
  ax_2 = fig.add_axes([x0 + width * left_percent + padding, y0, width * right_percent, height])

  # generate data to be smooth on graph
  comps = sorted(nx.connected_components(G), key=lambda x: len(x), reverse=True)
  comp = list(comps[0])
  data = np.zeros((G.order(), 1))
  for node in comp:
    data[node] = 0.6 + random.uniform(0.1, 0.3)
  node = random.choice(comp)
  data[node] = 0.95
  for comp in comps[1:]:
    for node in comp:
      data[node] = random.uniform(0.1, 0.3)

  row_range = np.arange(data.shape[0])
  ax_1.get_xaxis().set_visible(False)
  ax_1.invert_yaxis()
  ax_1.set_yticks(row_range + 0.5)
  ax_1.set_yticklabels(map(lambda x: str(x), row_range), fontsize=TICK_FONT_SIZE)
  if(column_index == 0):
    ax_1.set_ylabel("Node", fontsize=LABEL_FONT_SIZE)
  ax_1.set_title("Scores", fontsize=LABEL_FONT_SIZE)
  mappable = ax_1.pcolor(data, cmap=plt.cm.Reds, edgecolors='k', linewidth=LINE_WIDTH)

  ax_2.get_xaxis().set_visible(False)
  ax_2.get_yaxis().set_visible(False)
  plot_graph(G, ax_2, pos=pos)

  return data

# TODO this is a copy from man_reg_demo
def plot_graph(G, ax, pos=None):
  labels = {}
  for node in G.nodes():
    labels[node] = str(node)
  ax.get_xaxis().set_visible(False)
  ax.get_yaxis().set_visible(False)
  if pos is None:
    pos = nx.spring_layout(G)

  # draw scatter manually to add node borders
  xs = []
  ys = []
  for node_id, center in pos.items():
    xs.append(center[0])
    ys.append(center[1])
  nx.draw_networkx_edges(G, pos, ax=ax, width=LINE_WIDTH)
  # choose a light red for contrast with black
  rgba = plt.cm.Reds(0.4)
  ax.scatter(xs, ys, s=600, c=rgba, marker='o', edgecolor='black', linewidth=2.0, alpha=1.0)
  nx.draw_networkx_labels(G, pos, ax=ax, labels=labels, font_size=TICK_FONT_SIZE)

  return pos

def main():
  parser = argparse.ArgumentParser(description="Create a graphical abstract figure")
  parser.add_argument("out_png", type=str)
  args = parser.parse_args()

  data = np.array([0, 0.1, 0.8, 0.9, 0.9, 0.1])
  n_elem = data.shape[0]
  data = data.reshape((n_elem,1))

  fig_width = 1200
  fig_height = 900
  fig = plt.figure(figsize=(fig_width/DPI, fig_height/DPI), dpi=DPI)
  outer_ax = fig.gca()
  outer_coords = outer_ax.get_position().get_points().flatten()

  grid_row = 2
  grid_col = 3
  gs = mpl.gridspec.GridSpec(grid_row, grid_col,
    wspace=0.4, hspace=0.25)

  # create networks depicting pathways
  Gs = []
  for i in range(3):
    G = nx.generators.random_graphs.random_powerlaw_tree(5)
    edges = G.edges()
    edge = random.choice(edges)
    G.remove_edge(*edge)
    Gs.append(G)

  # visualize pairing of latent vectors with pathways -{{
  data = []
  for i in range(3):
    ax = plt.subplot(gs[1, i])
    ax.axis('off')
    datum = plot_latent_and_graph(Gs[i], fig, ax, i)
    data.append(datum)
  # }}-

  ax = plt.subplot(gs[0, :])
  aspect_ratio = fig_width / fig_height * grid_row
  n_latent = 4
  n_feature = 5
  n_data = len(data)
  for i in range(n_data, n_latent):
    data.append(None)
  plot_decomp(data, aspect_ratio, ax, n_latent=n_latent, n_feature=n_feature)

  left, bottom, right, top = outer_coords
  width = right - left
  height = top - bottom
  outer_ax = fig.add_axes([left, bottom, width, height])
  outer_ax.axis('off')
  outer_ax.annotate("", xytext=(0.7, 0.97), xy=(0.05, 0.50), arrowprops={
    'width': 5,
    'connectionstyle': 'angle3,angleA=0,angleB=90',
    'facecolor': 'grey'
  })
  outer_ax.annotate("", xytext=(0.7, 0.90), xy=(0.40, 0.50), arrowprops={
    'width': 5,
    'connectionstyle': 'angle3,angleA=0,angleB=90',
    'facecolor': 'grey'
  })
  outer_ax.annotate("", xytext=(0.7, 0.83), xy=(0.73, 0.50), arrowprops={
    'width': 5,
    'connectionstyle': 'angle3,angleA=90,angleB=0',
    'facecolor': 'grey'
  })

  fig.savefig(args.out_png, bbox_inches='tight')

if __name__ == "__main__":
  main()
