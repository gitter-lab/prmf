from . import *
import networkx as nx
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as colors

# https://www.infobyip.com/detectmonitordpi.php
DPI = 96
LINE_WIDTH = 2.0
LABEL_FONT_SIZE = 'xx-large'
TICK_FONT_SIZE = 'x-large'

def plot_vec(data, ax, aspect=3, yticklabels=None, vmin=None, vmax=None, log=False, column_index=0, n_columns=1):
  shape_tpl = data.shape
  data = data.reshape((shape_tpl[0], 1))
  if vmin is None:
    vmin = np.min(data)
  if vmax is None:
    vmax = np.max(data)

  # create labels centered at each heatmap block
  row_range = np.arange(shape_tpl[0])
  ax.get_xaxis().set_visible(False)
  ax.set_yticks(row_range + 0.5)
  if yticklabels is None:
    yticklabels = map(lambda x: str(x), row_range)
  ax.set_yticklabels(yticklabels, fontsize=TICK_FONT_SIZE)
  if column_index == 0:
    ax.set_ylabel("Node", fontsize=LABEL_FONT_SIZE)
  ax.set_title("Scores", fontsize=LABEL_FONT_SIZE)

  # only plot colorbar if column_index is the last column: n_columns - 1
  colorbar = None
  if log:
    mappable = ax.pcolor(data, cmap=plt.cm.Reds, edgecolors='k', linewidth=LINE_WIDTH, norm=colors.LogNorm(vmin=vmin, vmax=vmax), vmin=vmin, vmax=vmax)
    if column_index == n_columns - 1:
      colorbar = plt.colorbar(mappable, ax=ax, orientation='horizontal', ticks=[vmin, vmax], pad=0.05, fraction=0.08, aspect=aspect, format=mpl.ticker.LogFormatter(), use_gridspec=True)
      colorbar.set_ticks([vmin, vmax])
      colorbar.set_ticklabels(['{:.0e}'.format(vmin), '{:.0e}'.format(vmax)])
  else:
    mappable = ax.pcolor(data, cmap=plt.cm.Reds, edgecolors='k', linewidth=LINE_WIDTH, vmin=vmin, vmax=vmax)
    if column_index == n_columns - 1:
      colorbar = plt.colorbar(mappable, ax=ax, orientation='horizontal', ticks=[vmin, vmax], pad=0.05, fraction=0.08, aspect=aspect, use_gridspec=True)
  if column_index == n_columns - 1:
    colorbar.ax.tick_params(labelsize=TICK_FONT_SIZE)

def plot_graph(G, ax, pos=None, colormap=None, title='', title_fontsize=LABEL_FONT_SIZE, title_y=1.0, vmin=None, vmax=None, **kwargs):
  labels = {}
  for node in G.nodes():
    labels[node] = str(node)
  ax.get_xaxis().set_visible(False)
  ax.get_yaxis().set_visible(False)
  ax.set_title(title, fontsize=title_fontsize, y=title_y)
  if pos is None:
    pos = nx.spring_layout(G, k=1/np.sqrt(G.order()), iterations=200)

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

def plot_pathway_interactors(G, pathway_node, fig, ax, data, nodelist, title=None, title_fontsize=LABEL_FONT_SIZE, title_y=1.0, colormap=None, vmin=None, vmax=None, pos=None, **kwargs):
  if pathway_node not in G:
    raise FactorLibException("pathway_node must be a node in G")
  if colormap is None:
    colormap = plt.cm.Reds
  ax.get_xaxis().set_visible(False)
  ax.get_yaxis().set_visible(False)
  ax.axis('off')
  if vmin is None:
    vmin = np.min(data)
  if vmax is None:
    vmax = np.max(data)

  x0, y0, x1, y1 = ax.get_position().get_points().flatten()
  gs = mpl.gridspec.GridSpec(1, 2, width_ratios=[1,5], left=x0, bottom=y0, right=x1, top=y1)
  ax_1 = plt.Subplot(fig, gs[0])
  ax_2 = plt.Subplot(fig, gs[1])
  fig.add_subplot(ax_1)
  fig.add_subplot(ax_2)

  # - {{ ax_1 heatmap 
  # use node identifiers as heatmap labels but use node indices as graph labels because they dont fit well
  # associate node number with the identifier
  node_to_label = {}
  for i, node in enumerate(nodelist):
    node_to_label[node] = str(i)
  node_to_label[pathway_node] = pathway_node
  G = nx.relabel_nodes(G, node_to_label)
  label_list = []
  for i, node in enumerate(nodelist):
    label_list.append('{}: {}'.format(node, i))

  plot_vec(data, ax_1, yticklabels=label_list, vmin=vmin, vmax=vmax, log=True)
  # ax_1 heatmap }} -

  # - {{ ax_2 graph
  ax_2.get_xaxis().set_visible(False)
  ax_2.get_yaxis().set_visible(False)

  if pos is None:
    pos = nx.spring_layout(G, k=1/np.sqrt(G.order()), iterations=200)
  nx.draw_networkx_edges(G, pos, ax=ax_2, width=LINE_WIDTH)

  # draw scatter manually to add node borders
  xs = []
  ys = []
  for node_id, center in pos.items():
    if node_id != pathway_node:
      xs.append(center[0])
      ys.append(center[1])
  # choose a light red for contrast with black
  rgba = plt.cm.Reds(0.4)
  ax_2.scatter(xs, ys, s=600, c=rgba, marker='o', edgecolor='black', linewidth=LINE_WIDTH, alpha=1.0)

  # label all nodes other than pathway_node
  G.remove_node(pathway_node)
  pathway_node_center = pos.pop(pathway_node)
  nx.draw_networkx_labels(G, pos, ax=ax_2, font_size=TICK_FONT_SIZE)

  # label pathway_node
  ax_2.text(pathway_node_center[0], pathway_node_center[1], pathway_node, horizontalalignment='center', verticalalignment='center', fontsize=LABEL_FONT_SIZE, bbox=dict(facecolor=rgba, edgecolor='black', linewidth=LINE_WIDTH))

  ax_2.set_title(title, fontsize=title_fontsize, y=title_y)
  # ax_2 graph }} -

  return pos

def plot_latent_and_graph(G, fig, ax, data=None, nodelist=None, nodelabels=None,
  column_index=0, pos=None, colormap=None, max_col=30, log=True,
  title='', title_y=1.0, title_fontsize=LABEL_FONT_SIZE, 
  vmin=None, vmax=None, **kwargs):
  """
  Plot a vector and a graph
  """
  if data is not None and nodelist is None:
    raise FactorLibException("<data> and <nodelist> must be both set or both None")
  if colormap is None:
    colormap = plt.cm.Reds
  ax.get_xaxis().set_visible(False)
  ax.get_yaxis().set_visible(False)
  ax.axis('off')

  if data is None:
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
  else:
    # use provided data
    data = data.reshape((len(nodelist), 1))

  # split ax into 2+ subplots
  # if there are too many nodes to plot, separate heatmap into multiple pieces
  n_heatmaps = 1
  width_ratios = [1,5]
  if data.shape[0] > max_col:
    n_heatmaps = 2
    width_ratios = [1,1,1,4]
  n_col = len(width_ratios)

  x0, y0, x1, y1 = ax.get_position().get_points().flatten()
  gs = mpl.gridspec.GridSpec(1, n_col, width_ratios=width_ratios, left=x0, bottom=y0, right=x1, top=y1)
  axes = []
  for i in range(n_heatmaps):
    # skip every other ax in gridspec because it is used to make space for heatmap labels
    j = i*2
    ax = plt.Subplot(fig, gs[j])
    fig.add_subplot(ax)
    axes.append(ax)
  # add ax for graph
  ax = plt.Subplot(fig, gs[-1])
  fig.add_subplot(ax)
  axes.append(ax)

  # use node identifiers as heatmap labels but use node indices as graph labels because they dont fit well
  if nodelist is None:
    nodelist = sorted(G.nodes())
  node_to_ind = {}
  for i, node in enumerate(nodelist):
    node_to_ind[node] = i
  G = nx.relabel_nodes(G, node_to_ind)

  if vmin is None:
    vmin = np.min(data)
  if vmax is None:
    vmax = np.max(data)

  if nodelabels is None:
    nodelabels = []
    for i, node in enumerate(nodelist):
      nodelabels.append('{}: {}'.format(node, i))

  for i in range(n_heatmaps):
    data_start = 0 + max_col * i
    data_end = max_col * (i+1)
    data_this = data[data_start:data_end,0]
    yticklabels = nodelabels[data_start:data_end]
    plot_vec(data_this, axes[i], yticklabels=yticklabels, vmin=vmin, vmax=vmax, log=log, column_index=i, n_columns=n_heatmaps)

  ax = axes[-1]
  pos = plot_graph(G, ax, pos=pos, title=title, fontsize=title_fontsize, title_y=title_y)

  return pos

def split_title(title, n_lines=2):
  title_len = len(title)
  cur_len = 0
  new_words = []
  words = title.split()
  for word in words:
    word_len = len(word)
    new_words.append(word)
    cur_len += word_len
    if(cur_len) >= title_len / float(n_lines):
      new_words.append('\n')
  return ' '.join(new_words)
