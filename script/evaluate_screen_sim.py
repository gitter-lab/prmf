#!/usr/bin/env python
import sys, argparse
import os, os.path
import numpy as np
import sklearn.metrics
import factorlib as fl

def matching_id_to_ind(factor_id):
  return int(factor_id[1:])

def main():
  parser = argparse.ArgumentParser(description="""
Evaluate nmf_pathway.py results using simulated ground truth
""")
  parser.add_argument("--gene-by-latent", type=str, required=True, help=".csv file")
  parser.add_argument("--nodelist", type=argparse.FileType('r'), required=True)
  parser.add_argument("--true-seeds", type=argparse.FileType('r'), required=True)
  args = parser.parse_args()

  W_mat = np.genfromtxt(args.gene_by_latent, delimiter=",")
  nodelist = fl.parse_nodelist(args.nodelist)

  true_seed_fps = []
  for line in args.true_seeds:
    line = line.rstrip()
    true_seed_fps.append(line)

  true_seed_lists = []
  for true_seed_fp in true_seed_fps:
    seed_list = fl.parse_seedlist(true_seed_fp)
    true_seed_lists.append(seed_list)

  pathways_mat = fl.nodelists_to_mat(true_seed_lists, nodelist)
  
  # TODO probably error in normalize_num_pathways
  #W_mat, pathways_mat = fl.normalize_num_pathways(W_mat, pathways_mat)
  # TODO remove
  #np.savetxt("pathways_mat.csv", pathways_mat, delimiter=",")
  matching = fl.match(W_mat, pathways_mat)

  # 1) evaluate performance by a distance metric between latent factor and pathway
  # aucs are in range 0 to 1
  print("\t".join(["#" + "factor_id", "pathway_id", "auc"]))
  total_auc = 0
  n_aucs = 0
  for match in matching:
    factor_id, pathway_id, auc = match
    factor_id = matching_id_to_ind(factor_id)
    pathway_id = matching_id_to_ind(pathway_id)
    print("\t".join(map(str, match)))
    total_auc += auc
    n_aucs += 1
  print(total_auc / n_aucs)

if __name__ == "__main__":
  main()
