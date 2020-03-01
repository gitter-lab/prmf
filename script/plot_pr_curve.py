#!/usr/bin/env python
import argparse, sys
import os, os.path
import numpy as np
import sklearn
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import prmf

def matching_id_to_ind(factor_id):
  return int(factor_id[1:])

def main():
  parser = argparse.ArgumentParser(description="""
Evaluate NMF versus Pathway-Regularized Matrix Factorization by plotting PR curves on one figure.
""")
  parser.add_argument("--gene-by-latent-csvs", nargs="+", help=".csv files", required=True)
  parser.add_argument("--labels", nargs="+", help="parallel to --gene-by-latent-csvs", required=True)
  parser.add_argument("--nodelist", type=argparse.FileType('r'), required=True)
  parser.add_argument("--true-seeds", type=argparse.FileType('r'), required=True)
  parser.add_argument("--outdir", type=str, required=True)
  args = parser.parse_args()

  # parse inputs - {{
  W_mats = []
  colors = [] # TODO

  W_mats = list(map(lambda x: pd.read_csv(x, sep=",", header='infer', index_col=0).values, args.gene_by_latent_csvs))
  label_strs = list(map(lambda x: x + "; AUC={:0.3f}", args.labels))
  nodelist = prmf.parse_nodelist(args.nodelist)

  true_seed_fps = []
  for line in args.true_seeds:
    line = line.rstrip()
    true_seed_fps.append(line)

  true_seed_lists = []
  for true_seed_fp in true_seed_fps:
    seed_list = prmf.parse_seedlist(true_seed_fp)
    true_seed_lists.append(seed_list)

  pathways_mat = prmf.nodelists_to_mat(true_seed_lists, nodelist)
  # }} - parse inputs
  
  # reorganize <matching> so we can find each method's latent factor that best matches the ground truth
  pathway_to_latent_maps = []
  for i in range(len(W_mats)):
    matching = prmf.match(W_mats[i], pathways_mat)
    pathway_to_latent_map = {}
    for match in matching:
      factor_id_match, pathway_id_match, auc = match
      factor_id = matching_id_to_ind(factor_id_match)
      pathway_id = matching_id_to_ind(pathway_id_match)
      pathway_to_latent_map[pathway_id] = (factor_id, auc)
    pathway_to_latent_maps.append(pathway_to_latent_map)

  # plot Precision-Recall curves
  match_ind = 0
  method_to_avg_precision_vals = {}
  for i in range(len(W_mats)):
    method_to_avg_precision_vals[i] = []
  for pathway_id in range(pathways_mat.shape[1]):
    plt.clf()
    y_true = pathways_mat[:,pathway_id]
    true_fraction = np.sum(y_true) / y_true.shape[0]
    for i in range(len(pathway_to_latent_maps)):
      pathway_to_latent_map = pathway_to_latent_maps[i]
      factor_id, auc = pathway_to_latent_map[pathway_id]
      y_score = W_mats[i][:,factor_id]
      precision, recall, thresholds = sklearn.metrics.precision_recall_curve(y_true, y_score)
      method_to_avg_precision_vals[i].append(sklearn.metrics.average_precision_score(y_true, y_score))
      plt.plot(recall, precision, label=label_strs[i].format(auc), linewidth=2.0)
      #plt.plot(recall, precision, color=colors[i], label=label_strs[i].format(auc), linewidth=2.0)
    plt.plot(np.linspace(0,1,num=50), np.repeat(true_fraction, 50), label="Random; AUC={:0.3f}".format(true_fraction), linewidth=2.0)
    plt.xlabel('Recall', fontsize='x-large')
    plt.ylabel('Precision', fontsize='x-large')
    plt.ylim([0.0, 1.05])
    plt.xlim([0.0, 1.0])
    plt.title('Precision-Recall of PRMF and Friends', fontsize='xx-large')
    plt.legend()
    ofp = os.path.join(args.outdir, "fig{}.png".format(match_ind))
    plt.savefig(ofp, bbox_inches='tight')
    match_ind += 1

  # report the average of average precision for each method
  with open(os.path.join(args.outdir, 'avg_precision.txt'), 'w') as ofh:
    for i in range(len(W_mats)):
      avg_of_avg = np.mean(method_to_avg_precision_vals[i])
      label = args.labels[i]
      ofh.write("{}\t{}\n".format(label, avg_of_avg))

if __name__ == "__main__":
  main()
