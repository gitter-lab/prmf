#!/usr/bin/env python
import argparse
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
  ax_1.set_yticks(row_range + 0.5)
  ax_1.set_yticklabels(map(lambda x: str(x), row_range), fontsize=TICK_FONT_SIZE)
  if(column_index == 0):
    ax_1.set_ylabel("Node", fontsize=LABEL_FONT_SIZE)
  ax_1.set_title("Scores", fontsize=LABEL_FONT_SIZE)
  mappable = ax_1.pcolor(data, cmap=plt.cm.Reds, edgecolors='k', linewidth=LINE_WIDTH)

  ax_2.get_xaxis().set_visible(False)
  ax_2.get_yaxis().set_visible(False)
  plot_graph(G, ax_2, pos=pos)

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

  fig = plt.figure()
  fig = plt.figure(figsize=(900/DPI, 900/DPI), dpi=DPI)

  grid_row = 3
  grid_col = 3
  gs = mpl.gridspec.GridSpec(grid_row, grid_col)

  # create heatmap depicting input data -{{
  ax0 = plt.subplot(gs[0])
  data = np.random.rand(5,10)
  mappable = ax0.pcolor(data, cmap=plt.cm.Reds, linewidth=LINE_WIDTH)
  ax0.set_title("Features", fontsize=LABEL_FONT_SIZE)
  ax0.set_ylabel("Observations", fontsize=LABEL_FONT_SIZE)
  ax0.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
  ax0.tick_params(axis='y', which='both', left=False, right=False, labelleft=False)
  print(ax0.get_position().get_points())
  # }}-

  # create networks depicting pathways -{{
  ax1 = plt.subplot(gs[1])
  Gs = []
  for i in range(3):
    G = nx.generators.random_graphs.random_powerlaw_tree(5)
    edges = G.edges()
    edge = random.choice(edges)
    G.remove_edge(*edge)
    Gs.append(G)

  # stack a plot of a similar network
  # plots are stacked according to the order that the axes are created in
  xstep = 0.03
  ystep = 0.03
  bbox = ax1.get_position()
  x0, y0, x1, y1 = bbox.get_points().flatten()
  width = x1 - x0
  height = y1 - y0

  pos = [x0 + 2 * xstep, y0 + 2 * ystep, width, height]
  ax1_2 = fig.add_axes(pos)
  ax1_2.set_title("Pathways", fontsize=LABEL_FONT_SIZE)
  pos = [x0 + xstep, y0 + ystep, width, height]
  ax1_1 = fig.add_axes(pos)

  pos2 = plot_graph(Gs[2], ax1_2)
  pos1 = plot_graph(Gs[1], ax1_1)

  # redraw ax1 to put it on top
  plot_graph(Gs[0], ax1)
  pos = [x0, y0, width, height]
  ax1 = fig.add_axes(pos)
  pos0 = plot_graph(Gs[0], ax1)

  poses = [pos0, pos1, pos2]
  # }}-

  # visualize decomposition -{{
  # }}-

  # visualize pairing of latent vectors with pathways -{{
  for i in range(3):
    ax = plt.subplot(gs[2, i])
    ax.axis('off')
    plot_latent_and_graph(Gs[i], fig, ax, i, pos=poses[i])

  # }}-


  #G1_0
  #ax3_1 = inset_axes(parent_axes, width="30%", height=1, loc=3, border)

  fig.savefig(args.out_png)

if __name__ == "__main__":
  main()
