#!/usr/bin/env python
import sys, argparse
import os, os.path
import factorlib as fl
import numpy as np
import networkx as nx
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

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
    G = nx.read_graphml(pathway_fp)
    Gs.append(G)

  nodelist = fl.parse_nodelist(open(args.nodelist))

  latent_edge_diffs = fl.measure_smoothness(V, Gs, nodelist)

  for k in range(len(latent_edge_diffs)):
    csv_fh = open(os.path.join(args.outdir, 'latent{}.csv'.format(k)), 'w')
    hist_fp = os.path.join(args.outdir, 'latent{}.png'.format(k))

    edge_diffs = latent_edge_diffs[k]
    for edge_diff in edge_diffs:
      csv_fh.write("{}\n".format(edge_diff))

    # log parameter of plt.hist sets the y axis to log scale but not the x axis
    # change the x axis values to log scale
    edge_diffs = np.array(edge_diffs)

    plt.clf()
    n, bins, patches = plt.hist(edge_diffs)
    plt.xlabel('Edge differences')
    plt.ylabel('Frequency')
    plt.title('Smoothness of pathway edges\nmin = {:1.3f} ; max = {:1.3f}'.format(np.min(edge_diffs), np.max(edge_diffs)))
    plt.savefig(hist_fp)
