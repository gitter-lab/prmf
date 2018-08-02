#!/usr/bin/env python
import argparse
import os.path
import re
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
import factorlib as fl
from factorlib.plot import plot_latent_and_graph

# https://www.infobyip.com/detectmonitordpi.php
DPI = 96
LINE_WIDTH = 2.0
LABEL_FONT_SIZE = 'xx-large'
TICK_FONT_SIZE = 'x-large'

def truncate_graph(G, node_threshold=40):
  G = remove_isolated(G)
  conn_dict = {}
  edge_dict = {}
  comps = sorted(nx.connected_components(G), key=lambda x: len(x), reverse=True)
  largest_comp = comps[0]
  comp_size = len(largest_comp)
  while comp_size + len(conn_dict) > node_threshold:
    G = nx.subgraph(G, largest_comp)
    node_set = nx.minimum_node_cut(G)
    #if(len(node_set) != 1):
    #  raise FactorLibException("Graph is more highly connected than this function can handle")
    for node in node_set:
      edges = G.edges([node])
      G.remove_node(node)
      comps = sorted(nx.connected_components(G), key=lambda x: len(x), reverse=True)
      other_set = set()
      for comp in comps[1:]:
        for other_node in comp:
          other_set.add(other_node)
      conn_dict[node] = other_set
      largest_comp = comps[0]
      comp_size = len(largest_comp)
      edge_dict[node] = edges

  nodes = conn_dict.keys()
  nodes = sorted(nodes)
  for i, node in enumerate(nodes):
    group_id = "s{}".format(i)
    G.add_edge(node, group_id)
    edges = edge_dict[node]
    for edge in edges:
      other_node = None
      if edge[0] == node:
        other_node = edge[1]
      else:
        other_node = edge[0]
      if other_node in G:
        G.add_edge(node, other_node)

  return G, conn_dict, edge_dict

def remove_isolated(G):
  """
  Remove nodes with degree 0 from G
  """
  H = G.copy()
  node_to_degree = H.degree()
  for node, degree in node_to_degree.items():
    if degree == 0:
      H.remove_node(node)
  return H

def parse_mapping_file(fp):
  rv = {}
  with open(fp) as fh:
    for line in fh:
      words = line.split()
      rv[words[0]] = words[1]
  return rv

if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument("--pathway-mat")
  parser.add_argument("--pathway-obj")
  parser.add_argument("--nodelist")
  parser.add_argument("--truncate", default=False, type=bool, help="If True, truncate graph down to 50 nodes for visualization")
  parser.add_argument("--mapping-file", help="Node identifier mapping")
  parser.add_argument("--outdir")
  args = parser.parse_args()

  pathway_mat = np.genfromtxt(args.pathway_mat, delimiter=",")
  if(len(pathway_mat.shape) == 1):
    pathway_mat = pathway_mat.reshape(pathway_mat.shape[0], 1)
  latent_to_fp = fl.parse_pathway_obj(args.pathway_obj)
  latent_to_G = {}
  for k, fp in latent_to_fp.items():
    latent_to_G[k] = nx.read_graphml(fp)
  nodelist = fl.parse_nodelist(open(args.nodelist))
  mapping = parse_mapping_file(args.mapping_file)

  node_to_ind = {}
  for i, node in enumerate(nodelist):
    node_to_ind[node] = i

  for k, fp in latent_to_fp.items():
    G = latent_to_G[k]
    G = G.to_undirected()

    fig_width = 1200
    fig_height = 800
    fig = plt.figure(figsize=(fig_width/DPI, fig_height/DPI), dpi=DPI)

    fig.clf()
    ax = fig.gca()
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    ofp = os.path.join(args.outdir, "fig{}.png".format(k))

    # relabel nodes according to mapping and filter pathway nodes and their scores from latent factor
    vec = pathway_mat[:,k]
    node_to_score = {}
    vec_sub = []
    nodelist_sub = []
    for node in G.nodes():
      ind = node_to_ind[node]
      vec_sub.append(vec[ind])
      if node in mapping:
        nodelist_sub.append(mapping[node])
      else:
        nodelist_sub.append(node)
      node_to_score[node] = vec[ind]
    vec_sub = np.array(vec_sub)
    G = nx.relabel_nodes(G, mapping)
    for orig, new in mapping.items():
      if orig in node_to_score:
        # TODO assumes identifier lists orig and new do not overlap
        score = node_to_score.pop(orig)
        node_to_score[new] = score

    # TODO map fp to titles
    title_str = None
    bn = os.path.basename(os.path.splitext(fp)[0])
    bn = bn.split('_')[0]
    title_str = bn

    # DEBUG
    print(title_str)
    print('-' * len(title_str))
    print(vec_sub)
    print('')
    for n1, n2 in G.edges():
      s1 = node_to_score[n1]
      s2 = node_to_score[n2]
      diff = s1 - s2
      print(n1, n2, s1, s2, diff)
    print('')

    print("latent factor stats:")
    print("min: ", np.min(vec))
    print("max: ", np.max(vec))
    print("median: ", np.median(vec))
    print('')

    print("latent factor stats on pathway nodes:")
    print("min: ", np.min(vec_sub))
    print("max: ", np.max(vec_sub))
    print("median: ",np.median(vec_sub))
    print('')

    if args.truncate:
      G, conn_dict, edge_dict = truncate_graph(G)
      vec_sub, nodelist_sub = fl.filter_vec_by_graph(G, vec_sub, nodelist_sub)

    plot_latent_and_graph(G, fig, ax, data=vec_sub, nodelist=nodelist_sub, colormap=plt.cm.Reds, title=title_str, vmin=np.min(vec), vmax=np.max(vec))
    plt.savefig(ofp, bbox_inches='tight')
