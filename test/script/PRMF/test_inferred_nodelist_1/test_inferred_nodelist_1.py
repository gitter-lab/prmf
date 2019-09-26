#!/usr/bin/env python
import networkx as nx
import subprocess as sp
import numpy as np
import pandas as pd
import os
from scipy.stats import gamma

np.random.seed(seed=1)

# generate data from gamma distribution
m_samples = 100
n_genes = 1000
k_latent = 6

# X = U * V^T
# U is m_samples x k_latent
# V is n_genes x k_latent
shape_param = 5
U = gamma.rvs(shape_param, size=m_samples * k_latent).reshape(m_samples, k_latent)
V = gamma.rvs(shape_param, size=n_genes * k_latent).reshape(n_genes, k_latent)
X = U.dot(V.transpose())
nodelist = list(map(lambda x: "ENSP" + str(x), range(n_genes)))
X_df = pd.DataFrame(data=X, index=map(lambda x: str(x), range(m_samples)), columns=nodelist)

# write data
X_df.to_csv('data.tsv', sep='\t', header=False, index=False)

# write nodelist
with open('nodelist.txt', 'w') as fh:
  fh.write('\n'.join(nodelist))

# generate pathways
manifold_fps = []
for k in range(k_latent):
  inds = np.where(V[:,k] > np.percentile(V[:,k], 95))[0]
  G = nx.generators.path_graph(len(inds))

  # add name attribute to nodes
  for i in range(len(inds)):
    G.node[i]['name'] = "ENSP" + str(inds[i])

  manifold_fp = 'graph{}.graphml'.format(k)
  nx.write_graphml(G, manifold_fp)
  manifold_fps.append(manifold_fp)

# test with nodelist passed explicity
sp.check_call(args=['nmf_pathway.py', "--data", "data.tsv", "--manifolds"] + manifold_fps + ['--node-attribute', 'name', "--nodelist", 'nodelist.txt', "--outdir", os.curdir, '--delimiter', '\t', '--seed', "1"], stdout=open('test_inferred_nodelist_1.1.out', 'w'), stderr=open('test_inferred_nodelist_1.1.err', 'w'))
sp.check_call(args=['diff', '-q', 'test_inferred_nodelist_1_expected_obj.txt', 'obj.txt'])
print('test with nodelist passed')

# test with inferred nodelist but measurements on all genes
X_df.to_csv('data_header.tsv', sep='\t', header=nodelist, index=False)

sp.check_call(args=['nmf_pathway.py', "--data", "data_header.tsv", "--manifolds"] + manifold_fps + ['--node-attribute', 'name', "--outdir", os.curdir, '--delimiter', '\t', '--seed', "1"], stdout=open(os.devnull, 'w'), stderr=open(os.devnull, 'w'))
sp.check_call(args=['diff', '-q', 'test_inferred_nodelist_1_expected_obj.txt', 'obj.txt'])
print('test without nodelist passed')
