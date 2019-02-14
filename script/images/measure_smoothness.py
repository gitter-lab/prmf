#!/usr/bin/env python
import sys, argparse
import os, os.path
import itertools as it
import numpy as np
import networkx as nx
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import factorlib as fl

def measure_smoothness(V, Gs, nodelist, node_attr='name'):
  """
  Parameters
  ----------
  V: numpy.array
    result from nmf_pathway which is a gene x latent factor matrix

  Gs: list of networkx.Graph
    graphs with a node attribute called <node_attr>, 'name' by default, which is present in <nodelist>
    list is of length k_latent, matching the size of the array V
    
  nodelist: list of str
    mapping of row index in V to a gene symbol

  Returns
  -------
  TODO
  latent_data : list of list of float
    The outer list contains data for each latent factor.
    The inner list contains scores for each edge in the network associated with that latent factor.
    The edge scores are the absolute differences of node scores for that edge.
    This is a more detailed expansion of the data which is summarized by the manifold regularization metric.
  """
  n_genes, k_latent = V.shape

  node_to_index = {}
  for i, node in enumerate(nodelist):
    node_to_index[node] = i

  latent_data = []
  for k in range(k_latent):
    G = Gs[k]
    G = fl.relabel_nodes(G, node_attr)

    edge_diffs = []
    for u,v in G.edges_iter():
      u_index = node_to_index[u]
      v_index = node_to_index[v]
      edge_diff = abs(V[u_index,k] - V[v_index,k])
      edge_diffs.append(edge_diff)

    comp_diffs = []
    comps = sorted(nx.connected_components(G), key=lambda x: len(x), reverse=True)
    comp = comps[0]
    for u,v in it.combinations(comp, 2):
      u_index = node_to_index[u]
      v_index = node_to_index[v]
      diff = abs(V[u_index,k] - V[v_index,k])
      comp_diffs.append(diff)

    all_diffs = []
    for u,v in it.combinations(G.nodes(), 2):
      u_index = node_to_index[u]
      v_index = node_to_index[v]
      diff = abs(V[u_index,k] - V[v_index,k])
      all_diffs.append(diff)

    latent_datum = {
      'edge': edge_diffs,
      'comp': comp_diffs,
      'all': all_diffs
    }

    latent_data.append(latent_datum)
  return latent_data

def report_smoothness(measures, measure_type):
  for k in range(len(latent_edge_diffs)):
    csv_fh = open(os.path.join(args.outdir, 'latent{}_{}.csv'.format(k, measure_type)), 'w')
    hist_fp = os.path.join(args.outdir, 'latent{}_{}.png'.format(k, measure_type))

    edge_diffs = latent_edge_diffs[k][measure_type]
    for edge_diff in edge_diffs:
      csv_fh.write("{}\n".format(edge_diff))

    # log parameter of plt.hist sets the y axis to log scale but not the x axis
    # change the x axis values to log scale
    edge_diffs = np.array(edge_diffs)
    min_str = 'min = {:1.3f}'.format(np.min(edge_diffs)) if edge_diffs.shape[0] > 0 else 'min = NaN'
    max_str = 'max = {:1.3f}'.format(np.max(edge_diffs)) if edge_diffs.shape[0] > 0 else 'max = NaN'

    plt.clf()
    n, bins, patches = plt.hist(edge_diffs)
    plt.xlabel('{} differences'.format(measure_type.capitalize()))
    plt.ylabel('Frequency')
    plt.title('Smoothness of pathway edges\n{} ; {}'.format(min_str, max_str))
    plt.savefig(hist_fp)

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description="""
Measure the smoothness of latent factors for a converged result of nmf_pathway.
Plot the results in a histogram.
""")
  parser.add_argument("--indir", help="The output directory of nmf_pathway", required=True)
  parser.add_argument("--nodelist", help="Node to index association used by the nmf_pathway run", required=True)
  parser.add_argument("--outdir", help="Directory to place histogram and csv", required=True)
  args = parser.parse_args()

  V_fp = os.path.join(args.indir, 'V.csv')
  V = np.genfromtxt(V_fp, delimiter=",")
  n_genes, k_latent = V.shape

  obj_fp = os.path.join(args.indir, 'obj.txt')
  latent_to_pathway_fp = fl.parse_pathway_obj(obj_fp)
  Gs = []
  for k in range(k_latent):
    pathway_fp = latent_to_pathway_fp[k]
    G = nx.read_graphml(pathway_fp).to_undirected()
    Gs.append(G)

  nodelist = fl.parse_nodelist(open(args.nodelist))

  latent_edge_diffs = measure_smoothness(V, Gs, nodelist)
  report_smoothness(latent_edge_diffs, 'edge')
  report_smoothness(latent_edge_diffs, 'comp')
  report_smoothness(latent_edge_diffs, 'all')
