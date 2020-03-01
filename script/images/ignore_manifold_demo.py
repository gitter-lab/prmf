#!/usr/bin/env python
import argparse
import networkx as nx
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import rc
import prmf
import prmf.plot
rc('text', usetex=True)

DPI = prmf.plot.DPI
EPSILON = np.finfo(np.float32).eps

def project(w, G):
  inds = []
  for node in G:
    inds.append(node)
  inds = sorted(inds)
  return w[inds]

def manifold_penalty(w, G):
  L = nx.laplacian_matrix(G).todense()
  w = project(w, G)
  return L.dot(w).dot(w)

def ignore_penalty(w, G):
  # node identifiers are indexes
  rv = 0.0
  w = project(w, G)
  for i in range(w.shape[0]):
    rv += 1 / (w[i] + 1)
  return rv

if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument("--outfile")
  args = parser.parse_args()

  fig = plt.figure()
  fig = plt.figure(figsize=(1000/DPI, 500/DPI), dpi=DPI)

  gs = mpl.gridspec.GridSpec(1, 5, width_ratios=[1,5,1,1,5])
  ax1 = plt.subplot(gs[0,0])
  ax2 = plt.subplot(gs[0,1])
  #ax3 = plt.subplot(gs[0,2])
  ax4 = plt.subplot(gs[0,3])
  ax5 = plt.subplot(gs[0,4])

  # graph defined on 2-
  order = 7
  nodelist = list(range(order))
  G = nx.path_graph(order)
  G.remove_node(0)
  G.remove_node(1)
  #G.add_edge(3,6)

  # no penalty because smooth
  data1 = np.array([0.8, 0.9, 0.0, 0.0, 0.0, 0.0, 0.0])
  data1[data1 == 0] = EPSILON
  vmin = 0.0
  vmax = 1.0
  prmf.plot.plot_vec(data1, ax1, aspect=3, vmin=vmin, vmax=vmax)

  man_pen = manifold_penalty(data1, G)[0,0]
  title1 = "$w^T L w = {:1.3f}$".format(man_pen)

  # introduce penalty term for above data
  man_pen = manifold_penalty(data1, G)[0,0]
  ign_pen = ignore_penalty(data1, G)
  title2 = r'$w^T L w + \sum_L \frac{1}{w_L + 1} = ' + '{:1.3f}'.format(man_pen + ign_pen) + r'$'

  title = title1 + '\n' + title2
  pos = prmf.plot.plot_graph(G, ax2, title=title, title_fontsize=24, title_y=1.08)
  #prmf.plot.plot_graph(G, ax3, pos=pos, title=title)

  # penalty because not smooth on 2-
  data2 = np.array([0.0, 0.1, 0.8, 0.9, 0.8, 0.9, 0.7])
  data2[data2 == 0] = EPSILON
  prmf.plot.plot_vec(data2, ax4, aspect=3, vmin=vmin, vmax=vmax)

  man_pen = manifold_penalty(data2, G)[0,0]
  ign_pen = ignore_penalty(data2, G)
  title = r'$w^T L w + \sum_L \frac{1}{w_L + 1} = ' + '{:1.3f}'.format(man_pen + ign_pen) + '$'
  prmf.plot.plot_graph(G, ax5, pos=pos, title=title, title_fontsize=24, title_y=1.08)

  plt.savefig(args.outfile, bbox_inches='tight')
