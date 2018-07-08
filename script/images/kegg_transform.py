#!/usr/bin/env python
import argparse
import os, os.path
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import numpy as np
import networkx as nx
import factorlib.plot

DPI = 96
LINE_WIDTH = 2.0
LABEL_FONT_SIZE = 'xx-large'
TICK_FONT_SIZE = 'x-large'

# Shifted colorbar matplotlib
# https://stackoverflow.com/questions/20144529/shifted-colorbar-matplotlib/20146989#20146989
class MidpointNormalize(Normalize):
  def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
    self.midpoint = midpoint
    Normalize.__init__(self, vmin, vmax, clip)

  def __call__(self, value, clip=None):
    x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
    return np.ma.masked_array(np.interp(value, x, y))

# TODO this is a copy
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

if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument("outdir")
  args = parser.parse_args()

  G = nx.Graph()
  G.add_edge(0,1)
  G.add_edge(1,2)
  G.add_edge(0,3)
  G.add_edge(3,4)
  G.add_edge(5,6)
  G.add_edge(6,7)

  pos = {}
  pos[0] = (.1, .8)
  pos[1] = (.4, .8)
  pos[2] = (.7, .8)
  pos[3] = (.4, .5)
  pos[4] = (.7, .5)
  pos[5] = (.1, .2)
  pos[6] = (.4, .2)
  pos[7] = (.7, .2)

  fig = plt.figure(figsize=(391/DPI, 163/DPI), dpi=DPI)
  ax = fig.gca()
  ax.axis('off')
  plot_graph(G, ax, pos=pos)

  kegg_transform_fp = os.path.join(args.outdir, 'kegg_transform.png')
  fig.savefig(kegg_transform_fp, bbox_inches='tight')

  fig = plt.figure(figsize=(800/DPI, 800/DPI), dpi=DPI)
  ax = fig.gca()
  L = np.array(nx.laplacian_matrix(G).todense())

  vmin = np.min(L)
  vmax = np.max(L)

  # create labels centered at each heatmap block
  row_range = np.arange(L.shape[0])
  ax.set_xticks(row_range + 0.5)
  yticklabels = map(lambda x: str(x), row_range)
  ax.set_xticklabels(yticklabels, fontsize=factorlib.plot.TICK_FONT_SIZE)
  ax.set_yticks(row_range + 0.5)
  yticklabels = map(lambda x: str(x), row_range)
  ax.set_yticklabels(yticklabels, fontsize=factorlib.plot.TICK_FONT_SIZE)
  ax.set_ylabel("Node", fontsize=24)
  ax.set_xlabel("Node", fontsize=24)
  ax.xaxis.set_label_position('top') 
  ax.xaxis.tick_top()

  norm = MidpointNormalize(midpoint=0)
  mappable = ax.pcolor(L, cmap=plt.cm.RdGy, edgecolors='k', linewidth=factorlib.plot.LINE_WIDTH, vmin=vmin, vmax=vmax, norm=norm)
  colorbar = plt.colorbar(mappable, ax=ax, orientation='horizontal', ticks=[vmin, vmax], pad=0.05, fraction=0.08, aspect=6, use_gridspec=True)
  colorbar.ax.tick_params(labelsize=factorlib.plot.TICK_FONT_SIZE)

  laplacian_fp = os.path.join(args.outdir, 'kegg_laplacian.png')
  fig.savefig(laplacian_fp, bbox_inches='tight')
