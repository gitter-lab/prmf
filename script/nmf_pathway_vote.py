#!/usr/bin/env python
import sys, argparse
import os, os.path
import networkx as nx
from factorlib import script_utils
import pandas as pd
import numpy as np

def main():
  parser = argparse.ArgumentParser(description="""
Analyze a traditional NMF run's results by looking at each latent vector and measuring the
mass assigned to pathway nodes.
""")
  parser.add_argument("--gene-by-latent", help="Gene x Latent matrix from factorization with row and column names. Column names are e.g. \"LV1\"")
  parser.add_argument("--pathways", help="Graphml files where the node attribute \"name\" contains the gene identifier which is used as row names in --gene-by-latent", nargs="+")
  parser.add_argument("--outdir")
  args = parser.parse_args()
  script_utils.log_script(sys.argv)

  gene_by_latent_df = pd.read_csv(args.gene_by_latent, index_col=0)
  n_gene, k_latent = gene_by_latent_df.shape
  pathways = list(map(lambda x: nx.read_graphml(x), args.pathways))
  p_pathways = len(pathways)

  all_genes = set(gene_by_latent_df.index)
  latent_to_pathway_scores = np.zeros((k_latent, p_pathways))
  for k in range(k_latent):
    for p in range(p_pathways):
      pathway_members = list(map(lambda x: x[1]['name'], pathways[p].nodes(data=True)))
      #pathway_members_set = set(pathway_members)
      #non_pathway_members = all_genes - pathway_members_set

      # use avg pathway node mass as the latent-pathway matching score
      score = np.mean(gene_by_latent_df.loc[pathway_members, "LV{}".format(k)])
      latent_to_pathway_scores[k,p] = score

  # TODO change use of score to use maximum weighted matching
  pathway_to_argmax_arr = np.argmax(latent_to_pathway_scores, axis=1)
  print(pathway_to_argmax_arr)

if __name__ == "__main__":
  main()
