#!/usr/bin/env python
import sys, argparse
import os, os.path
import csv
import numpy as np
import networkx as nx
import pandas as pd
from scipy.stats import gamma
from scipy.stats import dirichlet
EPSILON = np.finfo(np.float32).eps

def main(args):
  X, U, V, pathways, other_pathways = run_simulation(args)
  save_simulation(outdir, X, U, V, pathways, other_pathways)

def run_simulation(args):
  p_other_pathways = 970
  pathway_order = 300
  k_latent = 30
  m_samples = 1000
  n_genes = 20000
  alpha = 5 # parameter for gamma distribution
  dirichlet_param = np.ones(k_latent) / k_latent

  np.random.seed(args.seed)

  pathways = []
  all_pathway_members = []
  for k in range(k_latent):
    pathway_members = np.random.randint(0, n_genes, size=pathway_order)
    # TODO raises NetworkXError
    # TODO the tree model is easy to work with but is not realistic because a gene's expression 
    # may be regulated by many proteins which would require a node to have multiple parents
    G = nx.random_powerlaw_tree(pathway_order, gamma=3, seed=seed, tries=100000)
    mapping = {}
    for i in range(len(pathway_members)):
      mapping[i] = pathway_members[i]
    H = nx.relabel_nodes(G, mapping)
    pathways.append(H)
    all_pathway_members.append(pathway_members)

  # generate other pathways that the data is not generated from
  other_pathways = []
  for p in range(p_other_pathways):
    pathway_members = np.random.randint(0, n_genes, size=pathway_order)
    G = nx.random_powerlaw_tree(pathway_order, gamma=3, seed=seed, tries=100000)
    mapping = {}
    for i in range(len(pathway_members)):
      mapping[i] = pathway_members[i]
    H = nx.relabel_nodes(G, mapping)
    other_pathways.append(H)

  # use the defined pathways to latent vectors
  # TODO if k_latent is different than p_pathways
  V_cols = []
  for k in range(k_latent):
    pathway = pathways[k]
    pathway_members = all_pathway_members[k]

    root = pathway_members[0]
    node_to_val = {}
    node_to_val[root] = gamma.rvs(alpha, size=1)[0]
    # NOTE tree assumptions in this loop
    for edge in nx.bfs_edges(pathway, root):
      parent_val = node_to_val[edge[0]]
      node_to_val[edge[1]] = gamma.rvs(parent_val + alpha, size=1)[0]

    V_col = gamma.rvs(alpha, size=n_genes)
    for node, val in node_to_val.items():
      V_col[node] = val
    V_col = V_col.reshape(n_genes, 1)

    V_cols.append(V_col)

  V = np.concatenate(V_cols, axis=1)

  # model mixtures of latent vectors with Dirichlet 
  U = dirichlet.rvs(dirichlet_param, size=m_samples)

  X = np.dot(U, V.transpose())

  return X, U, V, pathways, other_pathways

def save_simulation(outdir, X, U, V, Gs, Gs_other):
  m, n = X.shape
  m2, k = U.shape

  X = pd.DataFrame(X,
    index=map(lambda x: "sample{}".format(x+1), range(m)),
    columns=map(lambda x: "ENSG{}".format(x+1), range(n))
  )
  X.to_csv(os.path.join(outdir, "X.csv"), sep=",", quoting=csv.QUOTE_NONNUMERIC)

  U = pd.DataFrame(U, 
    index=map(lambda x: "sample{}".format(x+1), range(m)), 
    columns=map(lambda x: "LV{}".format(x+1), range(k))
  )
  U.to_csv(os.path.join(outdir, "U.csv"), sep=",", quoting=csv.QUOTE_NONNUMERIC)

  V = pd.DataFrame(V,
    index=map(lambda x: "ENSG{}".format(x+1), range(n)),
    columns=map(lambda x: "LV{}".format(x+1), range(k))
  )
  V.to_csv(os.path.join(outdir, "V.csv"), sep=",", quoting=csv.QUOTE_NONNUMERIC)

  pathway_files = []
  for i, G in enumerate(Gs):
    for node in G.nodes():
      G.node[node]['name'] = "ENSG{}".format(node+1)
    pathway_file = os.path.join(outdir, "pathway{}.graphml".format(i))
    pathway_files.append(pathway_file)
    nx.write_graphml(G, pathway_file)

  for i, G in enumerate(Gs_other):
    for node in G.nodes():
      G.node[node]['name'] = "ENSG{}".format(node+1)
    pathway_file = os.path.join(outdir, "pathway{}.graphml".format(i+len(Gs)+1))
    pathway_files.append(pathway_file)
    nx.write_graphml(G, pathway_file)

  with open(os.path.join(outdir, 'pathways_file.txt'), 'w') as fh:
    fh.write('\n'.join(pathway_files))

if __name__ == "__main__":
  parser = argparse.ArgumentParser(description="""
Generative simulation for matrix factorization.

For X approx UV^T where X is shape m_samples by n_genes, simulate sets of relationships among 
genes, called "pathways", which describe the dependency structure of gamma random variables 
defined on the genes.

Write files X.csv, U.csv, V.csv, pathway0.graphml, pathway1.graphml, ...pathwayN.graphml where 
pathway0.graphml, ... , pathway(K-1).graphml are the ground truth pathways. Write file 
pathways_file.txt which is a newline delimited file of the filepaths of pathway0.graphml, 
... , pathwayN.graphml.
""")
  parser.add_argument('--outdir', required=True, help="Output directory")
  args = parser.parse_args()
  main(args)
