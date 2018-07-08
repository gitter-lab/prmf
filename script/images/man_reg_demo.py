#!/usr/bin/env python
import argparse
import networkx as nx
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import rc
import factorlib.plot
rc('text', usetex=True)

# https://www.infobyip.com/detectmonitordpi.php
DPI = 96
LINE_WIDTH = 2.0
LABEL_FONT_SIZE = 'xx-large'
TICK_FONT_SIZE = 'x-large'
mpl.rcParams['lines.markeredgewidth'] = 5.0
#'scatter.marker'

def manifold_regularization(w, G):
  L = nx.laplacian_matrix(G)
  return L.dot(w).dot(w)

def main():
  parser = argparse.ArgumentParser(description="Create a heatmap from tabular CSV data")
  parser.add_argument("csv_file", type=str)
  parser.add_argument("out_png", type=str)
  args = parser.parse_args()

  data = np.genfromtxt(args.csv_file, delimiter=",")
  shape_tpl = data.shape
  if len(shape_tpl) == 1:
    data = data.reshape(data.shape[0],1)
  data_flat = data.flatten()

  fig = plt.figure()
  fig = plt.figure(figsize=(800/DPI, 480/DPI), dpi=DPI)

  gs = mpl.gridspec.GridSpec(2, 4, width_ratios=[1,3,0,3], height_ratios=[10, 1])
  ax0 = plt.subplot(gs[0,0])
  # indices swapped to put G2 on the left of G1
  ax1 = plt.subplot(gs[0,3])
  ax2 = plt.subplot(gs[0,1])

  # create labels centered at each heatmap block
  # TODO use fl.plot_vec
  row_range = np.arange(shape_tpl[0])
  ax0.get_xaxis().set_visible(False)
  ax0.set_yticks(row_range + 0.5)
  ax0.set_yticklabels(map(lambda x: str(x), row_range), fontsize=TICK_FONT_SIZE)
  ax0.set_ylabel("Node", fontsize=LABEL_FONT_SIZE)
  ax0.set_title("Scores", fontsize=LABEL_FONT_SIZE)
  mappable = ax0.pcolor(data, cmap=plt.cm.Reds, edgecolors='k', linewidth=LINE_WIDTH)
  colorbar = plt.colorbar(mappable, ax=ax0, orientation='horizontal', ticks=[0.0, 1.0], pad=0.05, fraction=0.08, aspect=3, use_gridspec=True)
  colorbar.ax.tick_params(labelsize=TICK_FONT_SIZE)
  
  # create networks illustrating manifold regularization values
  G1 = nx.path_graph(shape_tpl[0])
  G2 = G1.copy()
  G2.remove_edge(1,2)
  G2.remove_edge(4,5)

  # plot networks and manifold regularization value
  title1 = 'Penalty:\n$w^T L_2 w = {:1.3f}$'.format(manifold_regularization(data_flat, G1))
  pos = factorlib.plot.plot_graph(G1, ax1, title=title1, title_fontsize=24, title_y=1.08)
  title2 = 'Penalty:\n$w^T L_1 w = {:1.3f}$'.format(manifold_regularization(data_flat, G2))
  factorlib.plot.plot_graph(G2, ax2, pos=pos, title=title2, title_fontsize=24, title_y=1.08)

  ax3 = plt.subplot(gs[1,0])
  ax3.axis('off')
  ax3.text(0.5, 0.5, '$w$', horizontalalignment='center', verticalalignment='center', transform=ax3.transAxes, fontsize=LABEL_FONT_SIZE)

  ax4 = plt.subplot(gs[1,1])
  ax4.axis('off')
  ax4.text(0.5, 0.5, '$L_1$', horizontalalignment='center', verticalalignment='center', transform=ax4.transAxes, fontsize=LABEL_FONT_SIZE)

  ax5 = plt.subplot(gs[1,3])
  ax5.axis('off')
  ax5.text(0.5, 0.5, '$L_2$', horizontalalignment='center', verticalalignment='center', transform=ax5.transAxes, fontsize=LABEL_FONT_SIZE)

  fig.savefig(args.out_png, bbox_inches='tight')

if __name__ == "__main__":
  main()
