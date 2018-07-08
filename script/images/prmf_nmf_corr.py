#!/usr/bin/env python
import argparse
import os.path
import re
import numpy as np
from scipy.stats import pearsonr
import matplotlib.pyplot as plt

def parse_pathway_obj(fp):
  """
  Extract mapping of latent factor to pathway fp
  """
  rv = {}
  regexp = re.compile(r'(\d+)\s*->\s*(.+)')
  with open(fp) as fh:
    for line in fh:
      line = line.rstrip()
      match_data = regexp.match(line)
      if match_data is not None:
        rv[int(match_data.group(1))] = match_data.group(2)
  return rv

def match_col(pathway_mat, other_mat, col=0):
  m1,n1 = pathway_mat.shape
  m2,n2 = other_mat.shape
  corr_vec = np.zeros(n2,)
  for j in range(n2):
    corr = np.linalg.norm(pathway_mat[:,col] - other_mat[:,j]) / np.linalg.norm(pathway_mat[:,col])
    corr_vec[j] = corr
  
  # max along columns shows which other_mat column best matches a fixed pathway_mat column
  argmax = np.argmin(corr_vec)

  # https://stackoverflow.com/questions/20265229/rearrange-columns-of-numpy-2d-array
  # NOTE repeats OK
  #col_order = argmax_vec
  #i = np.argsort(col_order)
  #rv = other_mat[:,i]

  return other_mat[:,argmax], corr_vec[argmax]

if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument("--pathway-mat")
  parser.add_argument("--pathway-obj")
  parser.add_argument("--nmf-mats", nargs='+')
  parser.add_argument("--outdir")
  args = parser.parse_args()

  fig = plt.gcf()
  ax = fig.gca()

  pathway_mat = np.genfromtxt(args.pathway_mat, delimiter=",")
  latent_to_fp = parse_pathway_obj(args.pathway_obj)

  m,n = pathway_mat.shape
  corrs = []
  nmf_vecs = []
  for j in range(n):
    nmf_mat = np.genfromtxt(args.nmf_mats[0], delimiter=",")
    vec, corr = match_col(pathway_mat, nmf_mat, col=j)
    nmf_vecs.append(vec)
    corrs.append(corr)

    #def make_col_vec(vec):
    #  m = vec.shape[0]
    #  return vec.reshape(m,1)
    #arr = np.concatenate(list(map(make_col_vec, [pathway_mat[:,0]] + nmf_vecs)), axis=1)

  # plot correlation barchart, 1 column for each nmf_mat
  j = 0
  outfile = os.path.join(args.outdir, 'col_{}.png'.format(j))
  plt.clf()
  xs = np.arange(len(corrs))
  ys = np.array(corrs)
  plt.title('Correspondence between PRMF and NMF')
  plt.ylabel('Minimum norm of NMF column')
  plt.xlabel('PRMF column')
  plt.bar(xs, ys)
  plt.savefig(outfile, bbox_inches='tight')

  #LINE_WIDTH = 2.0
  #mappable = ax.pcolor(arr, cmap=plt.cm.Reds, linewidth=LINE_WIDTH)
  #ax.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
  #ax.tick_params(axis='y', which='both', left=False, right=False, labelleft=False)
  #plt.savefig(args.outfile, bbox_inches='tight')

