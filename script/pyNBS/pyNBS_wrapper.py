#!/usr/bin/env python
from pyNBS import data_import_tools as dit
from pyNBS import network_propagation as prop
from pyNBS import pyNBS_core as core
from pyNBS import pyNBS_single
from pyNBS import consensus_clustering as cc
from pyNBS import pyNBS_plotting as plot
import os
import time
import numpy as np
import argparse
import sys
import networkx as nx
import scipy.sparse as sp
import pandas as pd

def parse_nodelist(fh):
  """
  Return a list of node identifiers which maps the list index to the identifier
  """
  rv = []
  for line in fh:
    line = line.rstrip()
    words = line.split()
    for word in words:
      rv.append(word)
  return rv

def embed_ids(all_ids, ids):
  """
  Construct a binary vector of length len(<all_ids>) with a 1 in the positions
  associated with each id in <ids>
  """
  row_vec = np.zeros(len(all_ids))
  missing = []
  for id in ids:
    try:
      index = all_ids.index(id)
      row_vec[index] = 1
    except ValueError:
      missing.append(id)

  row_vec_sparse = sp.csc_matrix(row_vec)
  return (row_vec_sparse, missing)

def parse_gene_lists(nodelist, gene_list_fps):
  """
  Parse gene lists into a sparse scipy array
  
  Parameters
  ----------
  nodelist : list of str
    gene identifiers that define the assignment of array indexes to genes

  gene_lists : list of 

  Returns
  -------
  mat : sparse csc matrix
  """
  # parse gene lists
  gene_lists = []
  for gene_path in gene_list_fps:
    with open(gene_path) as fh:
      gene_lists.append(parse_nodelist(fh))

  # verify gene lists present in ppi_db
  def get_row_vec_for_gene_list(gene_list):
    row_vec, missing = embed_ids(nodelist, gene_list)
    sys.stderr.write("missing {}/{} node identifiers: {}\n".format(len(missing), len(gene_list), ", ".join(missing)))
    return row_vec
  row_vecs = map(get_row_vec_for_gene_list, gene_lists)

  mat = sp.vstack(row_vecs)
  return mat

if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument("--nodelist", required=True)
  parser.add_argument("--gene-lists", nargs='+', required=True)
  parser.add_argument("--network", required=True)
  parser.add_argument("--k-latent", required=True, type=int)
  parser.add_argument("--outdir", required=True)
  args = parser.parse_args()

  save_args = {'outdir': args.outdir, 'job_name': 'test'}

  nodelist = parse_nodelist(open(args.nodelist))
  arr_sp = parse_gene_lists(nodelist, args.gene_lists)
  arr_df = pd.DataFrame(arr_sp.todense())
  arr_df.columns = nodelist

  # graph regularizer
  network = nx.read_graphml(args.network)
  knnGlap = core.network_inf_KNN_glap(network)

  # diffusion
  alpha = 0.7
  network_nodes = network.nodes()
  network_I = pd.DataFrame(np.identity(len(network_nodes)), index=network_nodes, columns=network_nodes)
  kernel = prop.network_propagation(network, network_I, alpha=alpha, symmetric_norm=True)

  # Run pyNBS core steps
  H_df, W, H = pyNBS_single.NBS_single(arr_df, knnGlap, propNet=network, propNet_kernel=kernel, k=args.k_latent, **save_args)
  np.savetxt(os.path.join(args.outdir, "W.csv"), W, delimiter=",")
  np.savetxt(os.path.join(args.outdir, "H.csv"), H, delimiter=",")
