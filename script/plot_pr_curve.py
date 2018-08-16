#!/usr/bin/env python
import argparse, sys
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
  parser.add_argument("--nmf-gene-by-latent", type=str, help=".csv file")
  parser.add_argument("--prmf-gene-by-latent", type=str, help=".csv file")
  parser.add_argument("--plier-gene-by-latent", type=str, help=".csv file")
  parser.add_argument("--nodelist", type=argparse.FileType('r'), required=True)
  parser.add_argument("--true-seeds", type=argparse.FileType('r'), required=True)
  parser.add_argument("--outdir", type=str, required=True)
  args = parser.parse_args()

  if args.nmf_gene_by_latent is None and args.prmf_gene_by_latent is None and args.plier_gene_by_latent is None:
    sys.stderr.write("At least one of --nmf-gene-by-latent, --prmf-gene-by-latent, or --plier-gene-by-latent must be set\n")
    sys.exit(20)

  # parse inputs - {{
  W_mats = []
  label_strs = []
  colors = []
  if args.nmf_gene_by_latent is not None:
    W_mats.append(np.genfromtxt(args.nmf_gene_by_latent, delimiter=","))
    label_strs.append("NMF; AUC={:0.3f}")
    colors.append('grey')
  if args.prmf_gene_by_latent is not None:
    W_mats.append(np.genfromtxt(args.prmf_gene_by_latent, delimiter=","))
    label_strs.append("PRMF; AUC={:0.3f}")
    colors.append('red')
  if args.plier_gene_by_latent is not None:
    W_mats.append(np.genfromtxt(args.plier_gene_by_latent, delimiter=","))
    label_strs.append("PLIER; AUC={:0.3f}")
    colors.append('blue')

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
  # }} - parse inputs
  
  # reorganize <matching> so we can find each method's latent factor that best matches the ground truth
  pathway_to_latent_maps = []
  for i in range(len(W_mats)):
    matching = fl.match(W_mats[i], pathways_mat)
    pathway_to_latent_map = {}
    for match in matching:
      factor_id_match, pathway_id_match, auc = match
      factor_id = matching_id_to_ind(factor_id_match)
      pathway_id = matching_id_to_ind(pathway_id_match)
      pathway_to_latent_map[pathway_id] = (factor_id, auc)
    pathway_to_latent_maps.append(pathway_to_latent_map)

  # plot Precision-Recall curves
  match_ind = 0
  for pathway_id in range(pathways_mat.shape[1]):
    plt.clf()
    y_true = pathways_mat[:,pathway_id]
    for i in range(len(pathway_to_latent_maps)):
      pathway_to_latent_map = pathway_to_latent_maps[i]
      factor_id, auc = pathway_to_latent_map[pathway_id]
      y_score = W_mats[i][:,factor_id]
      precision, recall, thresholds = sklearn.metrics.precision_recall_curve(y_true, y_score)
      plt.step(recall, precision, color=colors[i], where='pre', label=label_strs[i].format(auc), linewidth=2.0)
    plt.xlabel('Recall', fontsize='x-large')
    plt.ylabel('Precision', fontsize='x-large')
    plt.ylim([0.0, 1.05])
    plt.xlim([0.0, 1.0])
    plt.title('Precision-Recall of PRMF and Friends', fontsize='xx-large')
    plt.legend()
    ofp = os.path.join(args.outdir, "fig{}.png".format(match_ind))
    plt.savefig(ofp, bbox_inches='tight')
    match_ind += 1

if __name__ == "__main__":
  main()
