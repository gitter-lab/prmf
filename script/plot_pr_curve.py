#!/usr/bin/env python
import argparse
import os, os.path
import numpy as np
import sklearn
import matplotlib.pyplot as plt
import factorlib as fl

def matching_id_to_ind(factor_id):
  return int(factor_id[1:])

def main():
  parser = argparse.ArgumentParser(description="""
Evaluate NMF versus Pathway-Regularized Matrix Factorization by plotting PR curves on one figure.
""")
  parser.add_argument("--nmf-gene-by-latent", type=str, required=True, help=".csv file")
  parser.add_argument("--prmf-gene-by-latent", type=str, required=True, help=".csv file")
  parser.add_argument("--nodelist", type=argparse.FileType('r'), required=True)
  parser.add_argument("--true-seeds", type=argparse.FileType('r'), required=True)
  parser.add_argument("--outdir", type=str, required=True)
  args = parser.parse_args()

  W_mat_nmf = np.genfromtxt(args.nmf_gene_by_latent, delimiter=",")
  W_mat_prmf = np.genfromtxt(args.prmf_gene_by_latent, delimiter=",")
  nodelist = fl.parse_nodelist(args.nodelist)

  true_seed_fps = []
  for line in args.true_seeds:
    line = line.rstrip()
    true_seed_fps.append(line)

  true_seed_lists = []
  for true_seed_fp in true_seed_fps:
    with open(true_seed_fp, 'r') as fh:
      seed_list = fl.parse_nodelist(fh)
      true_seed_lists.append(seed_list)

  pathways_mat = fl.nodelists_to_mat(true_seed_lists, nodelist)
  
  matching_nmf = fl.match(W_mat_nmf, pathways_mat)
  matching_prmf = fl.match(W_mat_prmf, pathways_mat)

  # transform matching_nmf to a dict keyed by pathway_id for comparison
  matching_nmf_dict = {}
  for match in matching_nmf:
    factor_id, pathway_id, auc = match
    matching_nmf_dict[pathway_id] = match

  # plot Precision-Recall curves
  match_ind = 0
  for match in matching_prmf:
    factor_id, pathway_id_match, auc_prmf = match
    factor_id_prmf = matching_id_to_ind(factor_id)
    pathway_id = matching_id_to_ind(pathway_id_match)

    factor_id, pathway_id_match2, auc_nmf = matching_nmf_dict[pathway_id_match]
    factor_id_nmf = matching_id_to_ind(factor_id)

    y_score_prmf = W_mat_prmf[:,factor_id_prmf]
    y_score_nmf = W_mat_nmf[:,factor_id_nmf]
    y_true = pathways_mat[:,pathway_id]

    precision_prmf, recall_prmf, thresholds_prmf = sklearn.metrics.precision_recall_curve(y_true, y_score_prmf)
    precision_nmf, recall_nmf, thresholds_nmf = sklearn.metrics.precision_recall_curve(y_true, y_score_nmf)
    #auc = sklearn.metrics.average_precision_score(y_true, y_score)

    plt.clf()
    plt.step(recall_prmf, precision_prmf, color='red', where='post', label="PRMF; AUC={:0.3f}".format(auc_prmf), linewidth=2.0)
    plt.step(recall_nmf, precision_nmf, color='grey', where='post', label="NMF; AUC={:0.3f}".format(auc_nmf), linewidth=2.0)
    #plt.fill_between(recall, precision, step='post', alpha=0.2, color='b')
    plt.xlabel('Recall', fontsize='x-large')
    plt.ylabel('Precision', fontsize='x-large')
    plt.ylim([0.0, 1.05])
    plt.xlim([0.0, 1.0])
    # TODO print area in legend
    plt.title('Precision-Recall of PRMF Against NMF'.format(auc), fontsize='xx-large')
    plt.legend()
    ofp = os.path.join(args.outdir, "fig{}.png".format(match_ind))
    plt.savefig(ofp, bbox_inches='tight')

    match_ind += 1

if __name__ == "__main__":
  main()
