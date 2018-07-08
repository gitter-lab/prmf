#!/usr/bin/env python
import argparse
import random
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import networkx as nx
import factorlib.plot as flp

if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument('--outfile')
  args = parser.parse_args()

  fig_width = 1200
  fig_height = 600
  fig = plt.figure(figsize=(fig_width/flp.DPI, fig_height/flp.DPI), dpi=flp.DPI)

  gs = mpl.gridspec.GridSpec(1, 5, width_ratios=[1, 3, 3, 1, 3])

  # plot heatmap on 15 nodes
  ax = plt.subplot(gs[0,0])
  ax.invert_yaxis()
  data = np.array([
    0.5, 0.6, 0.5, 0.7, 0.6,
    0.8, 0.7, 0.9, 0.5, 0.6,
    0.2, 0.1, 0.05, 0.2, 0.1
  ])
  flp.plot_vec(data, ax, vmin=0.0, vmax=1.0)

  inds = [1,2,4]
  for i, ind in enumerate(inds):
    mapping = {}
    for j, k in enumerate(range(i*5, (i+1)*5)):
      mapping[j] = k
    ax = plt.subplot(gs[0,ind])
    n_nodes = 5
    G = nx.random_powerlaw_tree(n_nodes)
    G = nx.relabel_nodes(G, mapping)
    flp.plot_graph(G, ax)

  ax = plt.subplot(gs[0,3])
  ax.axis('off')
  ax.axvline(x=0.5, linewidth=4, color='black')

  plt.savefig(args.outfile, bbox_inches='tight')
