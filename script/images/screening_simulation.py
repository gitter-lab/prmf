import argparse
import os, os.path
import numpy as np
import networkx as nx
import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy.sparse as sp
import factorlib as fl
import factorlib.plot as flp
EPSILON = np.finfo(np.float32).eps
DPI = flp.DPI

def make_G():
  """
  G and pos are meant to look like an excerpt of a pathway which contains the RIG-I-like Receptor Signaling Pathway
  """
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

  return G, pos

if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument("--outdir")
  args = parser.parse_args()

  n_nodes = 8
  G, pos = make_G()
  adj = nx.adjacency_matrix(G)

  vec = np.array([1.0, 0, 0, 0, 0, 0, 1.0, 0])
  vec[vec == 0] = EPSILON
  nodelist = list(map(lambda x: str(x), range(n_nodes)))
  nodelabels = nodelist

  fig_width = 800
  fig_height = 600
  fig = plt.figure(figsize=(fig_width/DPI, fig_height/DPI), dpi=DPI)
  gs = mpl.gridspec.GridSpec(1, 2, width_ratios=[1, 5])
  # diffusion in practice is done over a (transformed) global PPI network, plot a hairball instead
  H = nx.generators.random_graphs.binomial_graph(50, 0.2)
  flp.plot_vec(vec, plt.subplot(gs[0,0]), vmin=0.0, vmax=1.0, log=False)
  flp.plot_graph(H, plt.subplot(gs[0,1]))
  raw_ofp = os.path.join(args.outdir, 'raw.png')
  plt.savefig(raw_ofp, bbox_inches='tight')

  M = vec.reshape((1,n_nodes))
  M = sp.csc_matrix(M)
  vec_dif = fl.diffusion(M, adj, alpha=0.6)
  vec_dif = np.array(vec_dif.todense()).reshape((n_nodes,))
  print(vec_dif)

  # TODO plot vec, plot heatmap, draw arrow from diffused data into larger heatmap
  plt.clf()
  fig = plt.figure(figsize=(fig_width/DPI, fig_height/DPI), dpi=DPI)
  outer_ax = fig.gca()
  outer_coords = outer_ax.get_position().get_points().flatten()
  gs = mpl.gridspec.GridSpec(1, 3, width_ratios=[1,1,4])
  ax_1 = plt.subplot(gs[0,0])
  ax_2 = plt.subplot(gs[0,2])
  left, bottom, right, top = outer_coords
  width = right - left
  height = top - bottom
  outer_ax = fig.add_axes([left, bottom, width, height])
  outer_ax.axis('off')

  # plot vec
  flp.plot_vec(vec_dif, ax_1, vmin=0.0, vmax=1.0)

  # embed vec as part of larger array
  arr = np.random.rand(6, n_nodes)
  arr[4,:] = vec_dif

  # plot heatmap
  mappable = ax_2.pcolor(arr, cmap=plt.cm.Reds, linewidth=flp.LINE_WIDTH)
  ax_2.set_title("Genes", fontsize=flp.LABEL_FONT_SIZE)
  ax_2.set_ylabel("Observations", fontsize=flp.LABEL_FONT_SIZE)
  ax_2.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
  ax_2.tick_params(axis='y', which='both', left=False, right=False, labelleft=False)

  # plot annotation
  #  'connectionstyle': 'angle3,angleA=90,angleB=0',
  outer_ax.annotate("", xytext=(0.18, 0.73), xy=(0.4, 0.75), arrowprops={
    'width': 5,
    'facecolor': 'grey'
  })

  dif_ofp = os.path.join(args.outdir, 'diffuse.png')
  plt.savefig(dif_ofp, bbox_inches='tight')
