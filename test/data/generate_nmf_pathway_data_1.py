#!/usr/bin/env python
# Generate random data just to test the PRMF CLI
import numpy as np
import random
import os, os.path
import sys

random.seed(1)
np.random.seed(1)

this_fp = os.path.abspath(__file__)
data_dir = os.path.dirname(this_fp)
kegg_id_fp = os.path.join(data_dir, 'kegg_selection_ids.txt')
kegg_ids = []
with open(kegg_id_fp, 'r') as fh:
  for line in fh:
    line = line.rstrip()
    kegg_ids.append(line)

m = 100
n = 1000
k = 6

# X \approx UV^T
U = np.random.rand(m,k)
V_T = np.random.rand(k,n)
X = np.dot(U, V_T)
np.savetxt(os.path.join(data_dir, 'data.csv'), X, delimiter=",")

nodelist = random.sample(kegg_ids, n)
with open(os.path.join(data_dir, 'nodelist.txt'), 'w') as fh:
  fh.write('\n'.join(nodelist))
