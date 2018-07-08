#!/usr/bin/env python
import sys
import argparse
import os.path
import matplotlib.pyplot as plt
import numpy as np
import networkx as nx
import factorlib as fl
import factorlib.plot

def score_pathway_neighbors(ppi, pathway, nodelist, vec):
  """
  Find nodes connected (in <ppi>) to <pathway> but not members of <pathway> which have high scores in <vec>
  
  Parameters
  ----------
  ppi : nx.Graph

  pathway : nx.Graph

  nodelist : list of str

  vec : np.array

  Returns
  -------
  node_to_score : dict
    mapping of node identifier to score in vec
    includes nodes which are first neighbors of <pathway> as keys
  """
  node_to_score = {}
  pathway_set = set()
  for node in pathway.nodes():
    pathway_set.add(node)

  pathway_neighbors = set()
  for edge in ppi.edges(pathway.nodes()):
    n0, n1 = edge
    if n0 not in pathway_set:
      pathway_neighbors.add(n0)
    if n1 not in pathway_set:
      pathway_neighbors.add(n1)

  node_to_ind = {}
  for i, node in enumerate(nodelist):
    node_to_ind[node] = i

  for node in pathway_neighbors:
    ind = node_to_ind.get(node)
    if ind is None:
      node_to_score[node] = 0
      sys.stderr.write("Node \"{}\" not found in nodelist\n".format(node))
    else:
      node_to_score[node] = vec[ind]

  return node_to_score

def filter_pathway_neighbors(vec, node_to_score, max_nodes=10):
  rv = {}
  node_score_pairs = sorted(node_to_score.items(), key=lambda pair: pair[1], reverse=True)
  for node, score in node_score_pairs[:max_nodes]:
    rv[node] = score
  return rv

if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument('--nodelist', required=True)
  parser.add_argument('--gene-by-latent', required=True)
  parser.add_argument('--opt-outfile', required=True)
  parser.add_argument('--ppi-network', help="PPI and pathway union graph stored as graphml", required=True)
  parser.add_argument('--latent', default=None, help="If provided, only run script on this latent factor", type=int)
  parser.add_argument('--outdir', required=True)
  args = parser.parse_args()

  nodelist = fl.parse_nodelist(open(args.nodelist))
  ppi_network = nx.read_graphml(args.ppi_network)
  gene_by_latent = np.genfromtxt(args.gene_by_latent, delimiter=",")
  k_to_pathway_fp = fl.parse_pathway_obj(args.opt_outfile)
  if(args.latent is not None):
    k_to_pathway_fp = {args.latent: k_to_pathway_fp[args.latent]}

  ofp = os.path.join(args.outdir, 'pathway_extension.out')
  ofh = open(ofp, 'w')
  for k, fp in k_to_pathway_fp.items():
    pathway = nx.read_graphml(fp)
    vec = gene_by_latent[:,k]
    node_to_score = score_pathway_neighbors(ppi_network, pathway, nodelist, vec)
    node_to_score = filter_pathway_neighbors(vec, node_to_score)

    bn = os.path.basename(fp)
    bn, ext = os.path.splitext(bn)
    ofh.write(bn + '\n')
    for node, score in node_to_score.items():
      ofh.write(','.join([node, '{:.3f}'.format(score)]) + '\n')
    ofh.write('\n')

    pathway_node = bn.split("_")[0]
    pathway_name = pathway_node
    if pathway_name == "hsa04010":
      pathway_name = "MAPK Signaling Pathway"
    title = 'Top 10 {}-Interacting Genes'.format(pathway_name)
    title = factorlib.plot.split_title(title)

    plt.clf()
    fig = plt.gcf()
    ax = fig.gca()
    G_out = nx.subgraph(ppi_network, list(node_to_score.keys()))
    vec_sub, nodelist_sub = fl.filter_vec_by_graph(G_out, vec, nodelist)
    for node in node_to_score.keys():
      G_out.add_edge(node, pathway_node)
    factorlib.plot.plot_pathway_interactors(G_out, pathway_node, fig, ax, vec_sub, nodelist_sub, vmin=np.min(vec), vmax=np.max(vec), title=title, title_y=1.05)
    fig_out_fp = os.path.join(args.outdir, 'fig{}.png'.format(k))
    plt.savefig(fig_out_fp, bbox_inches='tight')
