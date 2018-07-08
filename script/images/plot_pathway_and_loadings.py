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
  parser.add_argument("--mapping-file", help="Node identifier mapping")
  parser.add_argument("--outdir")
  args = parser.parse_args()

  pathway_mat = np.genfromtxt(args.pathway_mat, delimiter=",")
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

    vec = pathway_mat[:,k]
    inds = []
    for node in G.nodes():
      ind = node_to_ind[node]
      inds.append(ind)
    # this may permute items in nodelist but the order of inds and vec_sub is consistent with nodelist_sub
    vec_sub = vec[inds]

    # relabel
    # TODO is order of nodes retained after mapping?
    G = nx.relabel_nodes(G, mapping)

    node_to_score = {}
    nodelist_sub = []
    for node in G.nodes():
      node_to_score[node] = vec[ind]
      nodelist_sub.append(node)

    # TODO map fp to titles
    title_str = None
    bn = os.path.basename(os.path.splitext(fp)[0])
    bn = bn.split('_')[0]
    if bn == "hsa04010":
      title_str = "MAPK Signaling Pathway"
    else:
      title_str = bn

    # DEBUG
    print(title_str)
    for n1, n2 in G.edges():
      s1 = node_to_score[n1]
      s2 = node_to_score[n2]
      print(n1, n2, s1, s2)
    print(inds)
    print(vec_sub)

    print("vec stats:")
    print(np.min(vec))
    print(np.max(vec))
    print(np.median(vec))

    print("vec_sub stats:")
    print(np.min(vec_sub))
    print(np.max(vec_sub))
    print(np.median(vec_sub))

    G, conn_dict, edge_dict = truncate_graph(G)
    vec_sub, nodelist_sub = fl.filter_vec_by_graph(G, vec_sub, nodelist_sub)

    plot_latent_and_graph(G, fig, ax, data=vec_sub, nodelist=nodelist_sub, colormap=plt.cm.Reds, title=title_str, vmin=np.min(vec), vmax=np.max(vec))
    plt.savefig(ofp, bbox_inches='tight')
