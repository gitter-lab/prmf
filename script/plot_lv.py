#!/usr/bin/env python
import argparse, sys
import pandas as pd
import matplotlib as mpl
import os, os.path
import networkx as nx
import numpy as np
mpl.use('Agg')
import matplotlib.pyplot as plt
from factorlib import script_utils
from factorlib import string_db
from factorlib import ensembl

def count_incident_nodes(G_ppi, G_pathway, nodes):
  """
  Count the number of connections among <nodes> to any of the nodes in <G_pathway>. 
  Do not include <nodes> in the count.

  Parameters
  ----------
  G_ppi : nx.Graph
    the network we are counting connections in, e.g. the STRINGdb protein-protein interaction network
    
  G_pathway : nx.Graph 
    the network defining which nodes we are counting connections to

  nodes : list of hashable
    the set of nodes that we are counting connections from 
  """
  G_sub = G_ppi.subgraph(list(G_pathway.nodes()) + list(nodes))
  incident_edges = []
  # TODO some of the edges may be connecting two nodes which are both mapped to the same ENSG (but it is unlikely)
  for node in nodes:
    incident_edges += G_sub.edges(node)
  incident_nodes = set()
  for edge in incident_edges:
    incident_nodes.add(edge[0])
    incident_nodes.add(edge[1])
    # incident_nodes may not include the node if there are no incident edges
  for node in nodes:
    if node in incident_nodes:
      incident_nodes.remove(node)
  count = len(incident_nodes)
  return count

def main():
  parser = argparse.ArgumentParser(description="""
Plot the top gene loadings for a latent vector
""")
  parser.add_argument("--gene-by-latent", "-V", required=True, help="Gene x Latent matrix from PRMF")
  parser.add_argument("--latent-index", "-l", required=True, help="Number indicating column in gene x latent matrix to plot", type=int)
  parser.add_argument("--n-genes", "-g", default=20, help="Number of genes to include in plot", type=int)
  parser.add_argument("--mapping-file", "-m", required=True, help="Ensembl mapping file associating HGNC to ENSG, ENST, ENSP in a .tsv file format")
  parser.add_argument("--pathway-file", "-p", required=True, help="Pathway proteins represented by ENSP identifiers, network in .graphml file format")
  parser.add_argument("--ppi-network", "-n", required=True, help="STRINGdb in ABC format with proteins represented by ENSP identifiers")
  parser.add_argument("--outdir")
  args = parser.parse_args()
  script_utils.log_script(sys.argv)

  df = pd.read_csv(args.gene_by_latent, header='infer', index_col=0)
  n_gene, k_latent = df.shape
  G_pathway = nx.read_graphml(args.pathway_file)
  G_ppi = string_db.parse_string_fh(open(args.ppi_network))
  latent_vec = df.iloc[:,args.latent_index]

  # note top scoring genes
  latent_vec = latent_vec.sort_values(ascending=False)
  top_scoring_genes = []
  top_scoring_values = []
  for gene, value in latent_vec.iloc[:args.n_genes].items():
    top_scoring_genes.append(gene)
    top_scoring_values.append(value)

  # note pathway gene scores
  pathway_vals = []
  pathway_ensgs = []
  ensp_to_ensg_map = ensembl.parse_mapping_tsv(args.mapping_file, key_index=3, value_index=1)
  for ensp in G_pathway.nodes():
    ensg_set = ensp_to_ensg_map.get(ensp)
    if ensg_set is not None:
      ensg = sorted(ensg_set)[0]
      pathway_vals.append(latent_vec.loc[ensg])
      pathway_ensgs.append(ensg)
  print(pathway_vals)

  # note connection of top scoring non-pathway genes to pathway genes with STRINGdb
  # STRINGdb is defined on ENSP
  print("Number of pathway proteins: {}".format(G_pathway.order()))
  ensg_to_ensp_map = ensembl.parse_mapping_tsv(args.mapping_file, key_index=1, value_index=3)
  incident_node_counts = []
  for gene in top_scoring_genes:
    count = 0
    protein_set = ensg_to_ensp_map.get(gene)
    if protein_set is None:
      print("Gene {} not found".format(gene))
    else:
      count = count_incident_nodes(G_ppi, G_pathway, protein_set)
    incident_node_counts.append(count)
  # note baseline connection to spliceosome
  # cache: 6.66704127503065
  # this takes 700 seconds to evaluate
  #count_sum = 0
  #for node in G_ppi.nodes():
  #  count_sum += count_incident_nodes(G_ppi, G_pathway, [node])
  #avg_spliceosome_connection = count_sum / G_ppi.order()
  #print(avg_spliceosome_connection)
  avg_spliceosome_connection = 6.66704127503065

  # map top scoring gene ensg to hgnc for plot x-axis labels
  ensg_to_hgnc = ensembl.parse_mapping_tsv(args.mapping_file, key_index=1, value_index=0)
  for i, ensg in enumerate(top_scoring_genes):
    hgnc_set = ensg_to_hgnc.get(ensg)
    if hgnc_set is not None:
      hgnc = sorted(hgnc_set)[0]
    else:
      hgnc = ensg
    top_scoring_genes[i] = hgnc

  # map pathway gene ensg to hgnc for plot x-axis labels
  pathway_hgnc_or_ensgs = []
  for i, ensg in enumerate(pathway_ensgs):
    hgnc_set = ensg_to_hgnc.get(ensg)
    if hgnc_set is not None:
      hgnc = sorted(hgnc_set)[0]
    else:
      hgnc = ensg
    pathway_hgnc_or_ensgs.append(hgnc)

  # generate plot of top scoring genes in the latent vector
  plt.clf()
  x_pos = np.arange(len(top_scoring_genes))
  plt.bar(x_pos, top_scoring_values, align='center')
  plt.xticks(x_pos, top_scoring_genes, rotation='vertical')
  plt.ylabel('LV{} loading'.format(args.latent_index))
  plt.title('Genes Representing LV{}'.format(args.latent_index))
  plt.tight_layout()
  plt.savefig(os.path.join(args.outdir, 'lv{}_genes.png'.format(args.latent_index)))

  # TODO pathway name is hard coded
  pathway_name = 'Spliceosome'

  # generate plot of pathway gene scores TODO
  hgnc_val_tpls = sorted(zip(pathway_hgnc_or_ensgs, pathway_vals), key=lambda x: x[1], reverse=True)
  pathway_hgnc_sorted = list(map(lambda x: x[0], hgnc_val_tpls))
  pathway_vals_sorted = list(map(lambda x: x[1], hgnc_val_tpls))
  plt.clf()
  x_pos = np.arange(args.n_genes)
  plt.bar(x_pos, pathway_vals_sorted[:args.n_genes], align='center')
  plt.xticks(x_pos, pathway_hgnc_sorted[:args.n_genes], rotation='vertical')
  plt.ylabel('LV{} loading'.format(args.latent_index))
  plt.title('Top {} Pathway Scores in LV{}'.format(pathway_name, args.latent_index))
  plt.tight_layout()
  plt.savefig(os.path.join(args.outdir, 'lv{}_pathway_genes.png'.format(args.latent_index)))

  # generate plot of number of pathway genes that each top scoring gene is connected to
  plt.clf()
  x_pos = np.arange(len(top_scoring_genes))
  plt.bar(x_pos, incident_node_counts, align='center')
  plt.plot(np.arange(-1, len(top_scoring_genes)+1), np.repeat(avg_spliceosome_connection, len(x_pos)+2), 'r--')
  plt.xticks(x_pos, top_scoring_genes, rotation='vertical')
  plt.ylabel('Number of pathway connections')
  plt.title('Connectedness of top-scoring LV{} genes to {}'.format(args.latent_index, pathway_name))
  plt.tight_layout()
  plt.savefig(os.path.join(args.outdir, 'lv{}_connections.png'.format(args.latent_index)))

if __name__ == "__main__":
  main()
