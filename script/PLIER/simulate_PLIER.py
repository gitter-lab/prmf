#!/usr/bin/env python
import argparse
import scipy.stats as st
import numpy as np
import networkx as nx
import random
import os, os.path

# return graphml fp
def sample_to_graphml(outdir, sample):
  G = nx.Graph()
  for j in sample:
    G.add_node(str(j))
  pathway_fp = os.path.join(outdir, "pathway{}.graphml".format(p))
  nx.write_graphml(G, pathway_fp)
  return pathway_fp

def sample_to_vec(n_genes, sample):
  rv = np.zeros((n_genes,))
  for i in sample:
    rv[i] = 1
  return rv

# || Y - ZB || + ||Z - CU || + ||B|| + ||U||_L1
# Y is gene x obs
# Z is gene x latent
# B is latent x obs
# C is gene x pathway
# U is pathway x latent
if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument("--outdir", "-o", required=True)
  parser.add_argument("--seed", "-s", default=1, type=int)
  args = parser.parse_args()
  random.seed(args.seed)
  np.random.seed(seed=args.seed)
  n_genes = 5000
  # "In this simulation we provide PLIER with 1000 pathways of which only 30 are correct"
  p_pathways = 30
  q_random_pathways = 970
  l_pathway_per_latent = 2
  k_latent = 30
  m_obs = 300

  # generate Z
  alpha = 5
  Z = st.gamma.rvs(alpha, size=n_genes*k_latent).reshape(n_genes, k_latent)

  # derive C
  C_true = np.zeros((n_genes, k_latent))
  pathway_fps = []
  for p in range(k_latent):
    # "threshold value on the percentage of genes which belong to a hypothetical prior knowledge geneset"
    threshold_value = random.uniform(1, 10)
    percentile_value = np.percentile(Z[:,p], 100 - threshold_value)
    inds = np.where(Z[:,p] > percentile_value)[0]
    for ind in inds:
      C_true[ind,p] = 1
    pathway_fps.append(sample_to_graphml(args.outdir, list(inds)))
  C_random = np.zeros((n_genes, q_random_pathways))
  pop = list(range(n_genes))
  for q in range(q_random_pathways):
    threshold_value = random.uniform(0.01, 0.1)
    pathway_size = int(threshold_value * n_genes)
    sample = random.sample(pop, pathway_size)
    pathway_fps.append(sample_to_graphml(args.outdir, sample))
    vec = sample_to_vec(n_genes, sample)
    C_random[:,q] = vec
  C = np.concatenate((C_true, C_random), axis=1)

  # force columns of B to sum to 1
  B = st.beta.rvs(1, 1, size=k_latent*m_obs).reshape(k_latent, m_obs)
  B_col_sums = np.sum(B, axis=0)
  for m in range(m_obs):
    B[:,m] = B[:,m] / B_col_sums[m]

  # construct Y
  E = st.norm.rvs(0, 1, size=n_genes*m_obs).reshape(n_genes, m_obs)
  Y = np.dot(Z, B) + np.abs(E)

  # save results
  np.savetxt(os.path.join(args.outdir, "Y.csv"), Y, delimiter=",")
  np.savetxt(os.path.join(args.outdir, "Z.csv"), Z, delimiter=",")
  np.savetxt(os.path.join(args.outdir, "B.csv"), B, delimiter=",")
  np.savetxt(os.path.join(args.outdir, "C.csv"), C, delimiter=",")
  with open(os.path.join(args.outdir, "nodelist.txt"), 'w') as fh:
    fh.write('\n'.join(map(str, range(1,n_genes+1))))
  with open(os.path.join(args.outdir, "pathways_file.txt"), 'w') as fh:
    fh.write('\n'.join(pathway_fps))
