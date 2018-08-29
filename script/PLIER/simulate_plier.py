import argparse
import numpy as np
import networkx as nx
import random
import os, os.path

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
  args = parser.parse_args()
  random.seed(1)
  n_genes = 10
  sample_size = 3
  p_pathways = 6
  k_latent = 4
  m_obs = 6

  # generate C
  C = np.zeros((n_genes, p_pathways))
  pop = list(range(n_genes))
  pathway_fps = []
  for p in range(p_pathways):
    sample = random.sample(pop, sample_size)
    G = nx.Graph()
    for j in sample:
      G.add_node(str(j))
    pathway_fp = os.path.join(args.outdir, "pathway{}.graphml".format(p))
    pathway_fps.append(pathway_fp)
    nx.write_graphml(G, pathway_fp)
    vec = sample_to_vec(n_genes, sample)
    C[:,p] = vec

  # generate a U that picks two pathways per latent factor
  U = np.zeros((p_pathways, k_latent))
  for k in range(k_latent):
    ps = random.sample(range(p_pathways), 2)
    for p in ps:
      U[p,k] = random.uniform(0.5, 1.5)

  # construct Z
  Z = np.dot(C,U)

  # construct a simple B such that each observation measures two latent factors
  # assume m_obs > k_latent
  B = np.zeros((k_latent, m_obs))
  for m in range(m_obs):
    k = m % k_latent
    k2 = (m + 1) % k_latent
    B[k,m] = random.uniform(0.25, 0.75)
    B[k2,m] = random.uniform(0.25, 0.75)

  # construct Y
  Y = np.dot(Z, B)

  # save results
  np.savetxt(os.path.join(args.outdir, "Y.csv"), Y, delimiter=",")
  np.savetxt(os.path.join(args.outdir, "Z.csv"), Z, delimiter=",")
  np.savetxt(os.path.join(args.outdir, "B.csv"), B, delimiter=",")
  np.savetxt(os.path.join(args.outdir, "C.csv"), C, delimiter=",")
  np.savetxt(os.path.join(args.outdir, "U.csv"), U, delimiter=",")
  with open(os.path.join(args.outdir, "nodelist.txt"), 'w') as fh:
    fh.write('\n'.join(map(str, range(n_genes))))
  with open(os.path.join(args.outdir, "pathways_file.txt"), 'w') as fh:
    fh.write('\n'.join(pathway_fps))
