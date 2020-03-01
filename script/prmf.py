#!/usr/bin/env python
import argparse, sys
from argparse import RawTextHelpFormatter
import numpy as np
import scipy.optimize
import scipy.sparse as sp
from scipy.stats import multinomial
from sklearn.preprocessing import quantile_transform
from sklearn.model_selection import train_test_split
from sklearn.model_selection import KFold
import networkx as nx
import prmf
from prmf import prmf_args
import copy
import datetime
import math
import os, os.path
import pandas as pd
import random
import csv
import datetime
np.seterr(divide='raise')
EPSILON = np.finfo(np.float32).eps

# TODO remove
OUTDIR = None
PATHWAY_TO_SUPPORT = None
LAPLACIANS = []
NORMALIZED_LAPLACIANS = []

def invert_list_map(ll):
  rv = {}
  for i, x in enumerate(ll):
    rv[x] = i
  return rv

def find_mins(V, Ls):
  """
  For each column vector v in V, find the L in Ls which minimizes v^T L v
  """
  n_feature, k_latent = V.shape
  rv = -1 * np.ones((k_latent))
  for k in range(k_latent):
    v = V[:,k]
    min_pen = np.inf
    min_ind = None
    for i, L in enumerate(Ls):
      # TODO change the order of dot products for faster computation
      man_penalty = (v.transpose().dot(L)).dot(v)
      if(man_penalty < min_pen):
        min_pen = man_penalty
        min_ind = i
    rv[k] = min_ind
  return rv

def normalize_laplacian(L, support):
  D_to_minus_half = sp.dia_matrix(L.shape)
  values = np.zeros(L.shape[0])
  for ind in support:
    v = L[ind, ind]
    if v != 0:
      values[ind] = v ** (-1/2)
  D_to_minus_half.setdiag(values)
  L_normal = D_to_minus_half.dot(L.dot(D_to_minus_half))
  return L_normal

def print_latent_to_pathway_data(latent_to_pathway_data):
  """
  TODO
  ----
  dump data to multiple latent x pathway csvs instead of printing to stdout

  Parameters
  ----------
  latent_to_pathway_data : dict<int, <int,double>>
    see restrict
  """
  print("-----------------------------")
  print("Latent to pathway match data:")
  print("-----------------------------")
  indent_level = 1
  indent = indent_level * 2 * " "
  for latent_id, pathway_data in latent_to_pathway_data.items():
    # sort by best scoring match descending
    print("Latent vector: {}".format(latent_id))
    pathway_data = sorted(pathway_data, key=lambda x: x[1], reverse=True)
    for pathway_id, score in pathway_data:
      print(indent + "{}\t{}".format(pathway_id, score))

def score_latent_pathway_match(pathway_to_support, Ls, latent_vec, pathway_id):
  """
  Parameters
  ----------
  pathway_to_support : dict<int, list<int>>
    mapping from pathway graph index to a list of nodes that are in the graph
    nodes are identified by integers which correspond to positions in <latent_vec>
    TODO aka lapl_to_feat_inds

  latent_vec : (n_genes, 1)-dimensional np.array

  pathway_id : 
    TODO aka pathway_ind

  """
  support = pathway_to_support[pathway_id]
  L = Ls[pathway_id]
  L_normal = normalize_laplacian(L, support)
  latent_vec_unit = latent_vec / np.linalg.norm(latent_vec)

  score_mass = np.sqrt(np.sum(np.power(latent_vec_unit[support], 2))) # in [0,1]
  score_manifold = 1 - L_normal.dot(latent_vec_unit).dot(latent_vec_unit) # in [0,1]
  score = score_mass + score_manifold
  return score

def score_latent_pathway_match_global(latent_vec, pathway_id):
  """
  See score_latent_pathway_match
  Use of global variables here is to test compute savings
  """
  latent_vec_unit = latent_vec / np.linalg.norm(latent_vec)
  support = PATHWAY_TO_SUPPORT[pathway_id]
  L_normal = NORMALIZED_LAPLACIANS[pathway_id]

  score_mass = np.sqrt(np.sum(np.power(latent_vec_unit[support], 2))) # in [0,1]
  score_manifold = 1 - L_normal.dot(latent_vec_unit).dot(latent_vec_unit) # in [0,1]
  score = score_mass + score_manifold
  return score

def restrict(V, Ls, latent_to_pathway_data, lapl_to_feat_inds):
  """
  For each column vector index pair (v,i) : v = V[:,i], find a subset of the Laplacians 
  Ls := latent_to_pathway_data[i] which have the property that the nodes in the graph associated
  with L_j are relatively overrepresented in the mass of v:
    
    TODO (relative change)
    \sum_j w[j : L_{jj} != 1] / sum_j w_j > #{j : L_{jj} != 1} / #{j}

  Parameters
  ----------
  V : np.array

  latent_to_pathway_data : dict<int, <int, double>>
    key is a column index for V
    value is a tuple of
      Laplacian index in Ls
      score for that Laplacian

  lapl_to_feat_inds : list of list
    outer list has length equivalent to p_pathways
    inner list contains the gene indexes where that pathway is defined (which genes from all 
      ~20,000 it is defined on)

  Returns
  -------
  rv : dict<int, int>
    updated latent_to_pathway_data
  """
  n_feature, k_latent = V.shape
  percentile = 19.9
  rv = {}
  for k in range(k_latent):
    pathway_data = latent_to_pathway_data[k]
    pathway_inds = list(map(lambda x: x[0], pathway_data))
    if(len(pathway_inds) > 1):
      v = V[:,k] 
      scores = np.zeros(len(pathway_inds))
      for i, L_ind in enumerate(pathway_inds):
        # scores[i] = score_latent_pathway_match(lapl_to_feat_inds, Ls, v, L_ind)
        scores[i] = score_latent_pathway_match_global(v, L_ind)

      score_inds = np.where(scores > np.percentile(scores,percentile))[0]

      pathway_ind_score_tuples = []
      if len(score_inds) == 0:
        # then the scores are probably uniform so that the top 80 (= 100 - 19.9) percent cannot be identified
        # randomly select 80 percent of the candidates instead
        n_candidates = len(pathway_inds)
        pathway_inds_arr = np.array(pathway_inds)
        sample_size = math.ceil(n_candidates * (1 - percentile) / 100)
        pathway_inds_rv = np.random.choice(pathway_inds_arr, size=sample_size, replace=False)
        pathway_scores_rv = list(map(lambda x: scores[x], pathway_inds_rv))
        pathway_ind_score_tuples = list(zip(pathway_inds_rv, pathway_scores_rv))
      else:
        for score_ind in score_inds:
          lapl_ind = pathway_inds[score_ind]
          score = scores[score_ind]
          ind_score_tpl = (lapl_ind, score)
          pathway_ind_score_tuples.append(ind_score_tpl)

      rv[k] = pathway_ind_score_tuples
    else:
      # otherwise converged to final Laplacian
      rv[k] = pathway_data
  return rv

def init_latent_to_pathway_data(k_latent, Ls):
  rv = {}
  for i in range(k_latent):
    rv[i] = list(zip(list(range(len(Ls))), [1] * len(Ls)))
  return rv

# TODO update this to add pathways to a set and count the size of the set
def count_distinct_pathways(latent_to_pathway_data):
  rv = np.inf
  for ind, lapls in latent_to_pathway_data.items():
    rv = min(rv, len(lapls))
  return rv

def force_distinct_lapls(V, Ls, latent_to_pathway_data, k_to_feat_inds, gamma, delta):
  """
  Finalize association between latent vector and pathway.

  For each candidate pathway that remains, evaluate the manifold and ignore penalties of the objective function.
  Construct a matching so that a latent vector is associated with the pathway which minimizes these penalty terms.

  Parameters
  -------
  latent_to_pathway_data : dict
    mapping from latent vector index to candidate pathway indexes
    NOTE mutated in place

  Returns
  -------
  latent_to_pathway_data : dict
    mapping from latent vector index to final pathway index
  """
  G = nx.Graph()
  for k, pathway_data in latent_to_pathway_data.items():
    for pathway_datum in pathway_data:
      lapl = pathway_datum[0]
      L = Ls[lapl]
      manifold_penalty = L.dot(V[:,k]).dot(V[:,k])
      ignore_penalty = 0
      for k2, nz_inds in k_to_feat_inds.items():
        ignore_penalty = np.sum(np.power(V[nz_inds,k2] + 1, -1))
      # TODO include mass for matching?
      denom = gamma * manifold_penalty + delta * ignore_penalty
      weight = None
      if denom == 0:
        # set to 0 because this situation only occurs when the latent vector has 0s on
        # the nodes that the Laplacian is constructed from
        weight = 0
      else:
        weight = 1 / denom
      G.add_edge("k{}".format(k), "l{}".format(lapl), {'weight': weight})
  mate = nx.max_weight_matching(G)
  for n1, n2 in mate.items():
    k_node = None
    l_node = None
    if n1[0] == 'k':
      k_node = int(n1[1:])
      l_node = int(n2[1:])
    else:
      k_node = int(n2[1:])
      l_node = int(n1[1:])
    # TODO update score value 2
    latent_to_pathway_data[k_node] = [(l_node, 2)]
  return latent_to_pathway_data

def map_k_to_lapls(k_to_lapl_ind, Ws, Ds, Ls, lapl_to_feat_inds):
  k_to_W = {}
  k_to_D = {}
  k_to_L = {}
  k_to_feat_inds = {}
  for k, lapl_ind in k_to_lapl_ind.items():
    k_to_W[k] = Ws[lapl_ind]
    k_to_D[k] = Ds[lapl_ind]
    k_to_L[k] = Ls[lapl_ind]
    k_to_feat_inds[k] = lapl_to_feat_inds[lapl_ind]
  return k_to_W, k_to_D, k_to_L, k_to_feat_inds

def pathway_to_vec(X, G, nodelist, rel_weight=5):
  """
  Construct a (n_gene x 1) vector which is smooth on the manifold given by G.
  Fill non-pathway genes with their average observation from X.

  Parameters
  ----------
  X : np.array
    m_obs x n_genes matrix

  G : nx.Graph
    undirected graph associated with a pathway

  nodelist : list of str
    mapping of array index to node identifier/gene symbol

  rel_weight : float
    the relative importance of the pathway nodes to the non-pathway nodes
  """
  n_genes = len(nodelist)
  v = np.zeros((n_genes,))
  node_to_index = {}
  for i, node in enumerate(nodelist):
    node_to_index[node] = i
  for node in G.nodes():
    v[node_to_index[node]] = 1
  pathway_ind = (v == 1)
  off_pathway_ind = np.invert(pathway_ind)
  v[off_pathway_ind] = np.mean(X.transpose()[off_pathway_ind], axis=1)
  off_pathway_mass = np.sum(v[off_pathway_ind])
  pathway_mass = off_pathway_mass * rel_weight
  v[pathway_ind] = pathway_mass / np.sum(pathway_ind)
  v = v.reshape((n_genes, 1))
  v = v / np.linalg.norm(v)
  return v, pathway_ind

def nmf_init_u(X, v):
  """
  Solve X = U V^T by non-negative least squares where U and V have shared dimension k_latent = 1
  and U,V >= 0 and V is initialized by <v>
  """
  m_obs, n_genes = X.shape
  u_prime, residual = scipy.optimize.nnls(X.transpose(), v.flatten())
  u = u_prime / np.linalg.norm(u_prime) ** 2
  u = u.reshape((m_obs,1))
  return u, residual

def nmf_init_v(X, u):
  """
  Solve X = U V^T where U is given and k_latent = 1

  TODO
  ----
  would like to impose equality constraint on variables corresponding to pathway members:
      min       || X v / ||v||_2^2  - u ||
     v >= 0
    subject to  v_i = v_j forall i,j in V(G)
  """
  m_obs, n_genes = X.shape
  v_prime, residual = scipy.optimize.nnls(X, u.flatten())
  v = v_prime / np.linalg.norm(v_prime) ** 2
  v = v.reshape(n_genes,1)
  return v, residual

def nmf_manifold_vec_obj(X, U, V, k_to_L, k_to_feat_inds, gamma=1, delta=1):
  obj_recon = np.linalg.norm(X - U.dot(V.transpose()))

  # TODO normal
  normal = True
  obj_manifold = 0.0
  obj_ign = 0.0
  if normal:
    for k, L in k_to_L.items():
      v_unit = V[:,k] / np.linalg.norm(V[:,k])

      # TODO could reorganize here to save recomputing normalized Laplacian
      support = k_to_feat_inds[k]
      L_normal = normalize_laplacian(L, support)
      obj_manifold += L_normal.dot(v_unit).dot(v_unit)

      obj_ign += np.sum(np.power(v_unit[support] + 1, -1))
  else:
    for k, L in k_to_L.items():
      obj_manifold += L.dot(V[:,k]).dot(V[:,k])
    for k, support in k_to_feat_inds.items():
      obj_ign += np.sum(np.power(V[support,k] + 1, -1))

  obj_fro = np.sum(np.multiply(U, U))
  obj_fro = obj_fro
  
  obj = obj_recon + gamma * obj_manifold + delta * obj_ign + obj_fro
  obj_data = {
    'recon': obj_recon,
    'manifold': obj_manifold,
    'ignore': obj_ign,
    'fro': obj_fro,
    'gamma': gamma,
    'delta': delta,
    'obj': obj
  }
  return obj_data

def nmf_manifold_vec_update(X, U, V, k_to_W, k_to_D, k_to_L, k_to_feat_inds, n_steps=10, gamma=1.0, delta=1.0, i=0, verbose=False, norm_X=None, tradeoff=0.5):
  """
  Perform <n_steps> update steps with a fixed Laplacian matrix for each latent factor

  Parameters
  ----------
  X : np.array
    data to factor

  U : np.array
    previous setting of U to update

  V : np.array
    previous setting of V to update

  k_to_W : dict
    mapping of latent factor to weighted adjacency matrix

  k_to_D : dict
    mapping of latent factor to diagonal matrix that is the sum of W along a row (or column)

  k_to_L : dict
    mapping of latent factor to L = D - W

  n_steps : int
    number of update steps to perform

  gamma : float
    relative importance of manifold regularization term

  delta : float
    relative importance of ignoring manifold penalty

  i : int
    number of previous iterations

  verbose : bool
    if True, print objective function value after each iteration

  norm_X : float or None
    stored value of the norm of X
  """
  obj_data = None
  m, k_latent = U.shape
  n, k_latent = V.shape
  for n_step in range(n_steps):
    U_up_num = X.dot(V)
    U_up_denom = U.dot((V.transpose().dot(V))) + U
    U = np.multiply(U, np.divide(U_up_num, U_up_denom, out=np.ones_like(U_up_num), where=U_up_denom!=0)) # 0 / 0 := 1

    V_up_num_recon = X.transpose().dot(U)
    V_up_denom_recon = V.dot((U.transpose().dot(U)))

    # update each column vector of V separately to accomodate different Laplacians
    V_up_num_man = np.zeros((n, k_latent))
    V_up_denom_man = np.zeros((n, k_latent))
    V_up_num_ign = np.zeros((n, k_latent))
    for k in range(k_latent):
      W = k_to_W[k]
      D = k_to_D[k]
      V_up_num_man[:,k] = gamma * W.dot(V[:,k])
      V_up_denom_man[:,k] = gamma * D.dot(V[:,k])

      nz_inds = k_to_feat_inds[k]
      V_up_num_ign[nz_inds,k] = delta * np.power(V[nz_inds,k] + 1, -2)

    V_up_num = V_up_num_recon + (V_up_num_man + V_up_num_ign)
    V_up_denom = V_up_denom_recon + V_up_denom_man
    V_up_denom[V_up_denom < EPSILON] = EPSILON
    V = np.multiply(V, np.divide(V_up_num, V_up_denom, out=np.ones_like(V_up_num), where=V_up_denom!=0))
    V[V < EPSILON] = EPSILON

    obj_data = nmf_manifold_vec_obj(X, U, V, k_to_L, k_to_feat_inds, gamma=gamma, delta=delta)
    print(i+n_step+1, obj_data['obj'])
    if(verbose):
      print(obj_data)

  return U, V, obj_data

def nmf_manifold_vec_update_normal(X, U, V, k_to_W, k_to_D, k_to_L, k_to_feat_inds, n_steps=10, gamma=1.0, delta=1.0, i=0, verbose=False, norm_X=None, tradeoff=0.5):
  """
  See nmf_manifold_vec_update ; this uses the normalized Laplacian instead
  """
  obj_data = None
  m, k_latent = U.shape
  n, k_latent = V.shape
  for n_step in range(n_steps):
    U_up_num = X.dot(V)
    U_up_denom = U.dot((V.transpose().dot(V))) + U
    U = np.multiply(U, np.divide(U_up_num, U_up_denom, out=np.ones_like(U_up_num), where=U_up_denom!=0)) # 0 / 0 := 1

    V_up_num_recon = X.transpose().dot(U)
    V_up_denom_recon = V.dot((U.transpose().dot(U)))

    # update each column vector of V separately to accomodate different Laplacians
    # TODO need to use v_unit in manifold term
    V_up_num_man = np.zeros((n, k_latent))
    V_up_denom_man = np.zeros((n, k_latent))
    V_up_num_ign = np.zeros((n, k_latent))
    for k in range(k_latent):
      W = k_to_W[k]
      D = k_to_D[k]
      D_to_minus_half = D.power(-1/2)
      v_norm_sq_inv = 1/(np.linalg.norm(V[:,k])**2)
      V_up_num_man[:,k] = gamma * v_norm_sq_inv * D_to_minus_half.dot(W.dot(D_to_minus_half)).dot(V[:,k])
      V_up_denom_man[:,k] = gamma * v_norm_sq_inv * V[:,k]

      nz_inds = k_to_feat_inds[k]
      V_up_num_ign[nz_inds,k] = delta * v_norm_sq_inv * np.power(V[nz_inds,k] + 1, -2)

    V_up_num = V_up_num_recon + (V_up_num_man + V_up_num_ign)
    V_up_denom = V_up_denom_recon + V_up_denom_man
    V_up_denom[V_up_denom < EPSILON] = EPSILON
    V = np.multiply(V, np.divide(V_up_num, V_up_denom, out=np.ones_like(V_up_num), where=V_up_denom!=0))
    V[V < EPSILON] = EPSILON

    obj_data = nmf_manifold_vec_obj(X, U, V, k_to_L, k_to_feat_inds, gamma=gamma, delta=delta)
    print(i+n_step+1, obj_data['obj'])
    if(verbose):
      print(obj_data)

  return U, V, obj_data

def nmf_manifold_vec_update_tradeoff(X, U, V, k_to_W, k_to_D, k_to_L, k_to_feat_inds, n_steps=10, i=0, verbose=False, norm_X=None, tradeoff=0.5, gamma=1, delta=1):
  """
  See nmf_manifold_vec_update; this version sets gamma and delta for the _next_ gradient descent step 
  so that delta = gamma = (1 - tradeoff) * obj_recon / (tradeoff * obj_manifold)

  Parameters
  ----------
  tradeoff : float
    value in [0,1] representing relative importance of reconstruction error to manifold regularization penalty.
    alternative to gamma and delta. 1 means only use reconstruction error.
  """
  obj_data = None
  m, k_latent = U.shape
  n, k_latent = V.shape
  for n_step in range(n_steps):
    U_up_num = X.dot(V)
    U_up_denom = U.dot((V.transpose().dot(V))) + U
    U = np.multiply(U, np.divide(U_up_num, U_up_denom, out=np.ones_like(U_up_num), where=U_up_denom!=0)) # 0 / 0 := 1

    V_up_num_recon = X.transpose().dot(U)
    V_up_denom_recon = V.dot((U.transpose().dot(U)))

    # update each column vector of V separately to accomodate different Laplacians
    V_up_num_man = np.zeros((n, k_latent))
    V_up_denom_man = np.zeros((n, k_latent))
    V_up_num_ign = np.zeros((n, k_latent))
    for k in range(k_latent):
      W = k_to_W[k]
      D = k_to_D[k]
      V_up_num_man[:,k] = gamma * W.dot(V[:,k])
      V_up_denom_man[:,k] = gamma * D.dot(V[:,k])

      nz_inds = k_to_feat_inds[k]
      V_up_num_ign[nz_inds,k] = delta * np.power(V[nz_inds,k] + 1, -2)

    V_up_num = V_up_num_recon + (V_up_num_man + V_up_num_ign)
    V_up_denom = V_up_denom_recon + V_up_denom_man
    V_up_denom[V_up_denom < EPSILON] = EPSILON
    V = np.multiply(V, np.divide(V_up_num, V_up_denom, out=np.ones_like(V_up_num), where=V_up_denom!=0))
    V[V < EPSILON] = EPSILON

    obj_data = nmf_manifold_vec_obj(X, U, V, k_to_L, k_to_feat_inds, gamma=gamma, delta=delta)
    # update gamma and delta such that at the next iteration the reconstruction error contributes 
    # <tradeoff> portion of the objective function and the manifold regularization contributes 
    # 1 - <tradeoff> portion
    denom = (tradeoff * obj_data['manifold'])
    if denom == 0:
      # then unscaled manifold penalty is near 0 anyway, define x / 0 := 1
      gamma = 1
    else:
      gamma = ((1 - tradeoff) * obj_data['recon']) / denom
    delta = gamma

    print(i+n_step+1, obj_data['obj'])
    if(verbose):
      print(obj_data)

  return U, V, obj_data, gamma, delta

def nmf_pathway(X, Gs, gamma=1.0, delta=1.0, tradeoff=None, k_latent=6, tol=1e-3, max_iter=1000, nodelist=None, modulus=10, U_init=None, V_init=None, verbose=False):
  """
  Solve an optimization problem of the form
    min ||X - UV^T|| + 
      gamma * sum_k min_i V[:,k]^T Ls[i] V[:,k] + 
      delta * sum_k sum_{i | i in G_k} 1 / V[i,k] + 
      ||U||_F^2

  where Ls[i] is the Laplacian matrix associated with Gs[i],
  G_k is the manifold associated with latent factor k
  X has shape (n_obs, n_features),
  U has shape (n_obs, n_latent),
  V has shape (n_feature, n_latent)
  
  Parameters
  ----------
  X : np.array
    matrix of shape (m_obs, n_features)

  Gs : list of nx.Graph
    set of sparse manifolds where the data matrix X is assumed to be smooth with respect to 
    one or more manifolds present in this set. Graphs may set the 'weight' attribute in [0,1] to
    indicate confidence in the presence of the edge. 1 is high.

  gamma : float
    trade-off between reconstruction error and manifold regularization penalty

  delta : float
    regularization parameter for penalty for ignoring manifold

  tradeoff : None, float
    if not None, triggers automatic setting of gamma and delta after each iteration so that 
    the relative importance of the reconstruction error and the manifold regularization penalties
    (including penalty for ignoring the manifold) is at a fixed proportion:
      tradeoff * reconstruction_error + (1 - tradeoff) * manifold_regularization_penalties

  k_latent : int
    number of latent factors to decompose X into: for X \\approx UV^T, U is shape (m_obs, k_latent)
    and V^T is shape (k_latent, n_features)

  tol : float
    objective function convergence tolerance. iterative method has converged when successive
    solutions do not improve by more than <tol>

  max_iter : int
    maximum number of iterations. function terminates and returns the current solution 
    after this number of iterations.

  nodelist : list of str
    node identifiers for nodes in any graph in G. fixes an order of the nodes so that the graphs in G
    may be translated to Laplacian matrices where a node in the graph is mapped to a particular index
    according to its position in the nodelist.

  modulus : int
    number of iterations to retain active Laplacian and candidate Laplacians

  U_init : np.array
    array with shape (n_obs, n_latent) used to initialize PRMF

  V_init : np.array
    array with shape (n_feature, n_latent) used to initialize PRMF

  verbose : bool
    if true, increase printing to stdout

  Returns
  -------
  U : np.array

  V : np.array

  obj_data : dict
    information about the objective function at termination has keys
    'recon': reconstruction error : ||X - WH||
    'manifold': manifold regularization penalty: gamma * sum_k min_i W[:,k]^T Ls[i] W[:,k]
    'obj': objective function value: sum of 'recon' and 'manifold'
  """
  global PATHWAY_TO_SUPPORT
  global LAPLACIANS
  global NORMALIZED_LAPLACIANS
  # rescale to prevent underflow errors
  #alpha = 1 / np.min(X[X != 0])
  alpha = 1
  X = alpha * X
  norm_X = np.linalg.norm(X)
  print('norm(X) = {}'.format(norm_X))

  # TODO note gamma and delta rescaling
  gamma = gamma * norm_X / k_latent
  delta = delta * 10 / norm_X

  m,n = X.shape
  # 1 - [0,1) \in (0,1] ; need strictly positive
  # TODO other initialization strategies
  if U_init is None:
    U_init = 3 * (1 - np.random.rand(m,k_latent))
  if V_init is None:
    V_init = 3 * (1 - np.random.rand(n,k_latent))
  U = U_init
  V = V_init
  if U.shape != (m, k_latent):
    raise ValueError('Invalid U_init with shape {} != (m_obs, k_latent) = {}'.format(U_init.shape, (m, k_latent)))
  if V.shape != (n, k_latent):
    raise ValueError('Invalid V_init with shape {} != (n_feature, k_latent) = {}'.format(V_init.shape, (n, k_latent)))

  # TODO nodelist cant be None, initialization?
  node_to_index = {}
  for i, node in enumerate(nodelist):
    node_to_index[node] = i

  # preprocess networks and nodelist
  # the nodelist specifies the entire universe of identifiers for this run
  # remove nodes (and incident edges) which are not present in the nodelist
  # TODO warn? move this description to --help
  for i, G in enumerate(Gs):
    Gs[i] = G.subgraph(nodelist)

  Ws = []
  Ds = []
  Ls = []
  i = 0
  lapl_to_feat_inds = {}
  for i, G in enumerate(Gs):
    W = nx.adjacency_matrix(G, nodelist=nodelist)
    data = W.sum(axis=0)
    offsets = np.array([0])
    D = sp.dok_matrix(sp.dia_matrix((data, offsets), shape=(n,n)))
    L = D - W
    Ws.append(W)
    Ds.append(D)
    Ls.append(L)

    # features in the Laplacian matrix that have non-zero entries on the diagonal
    feat_inds = list(map(lambda x: node_to_index[x], G.nodes()))
    lapl_to_feat_inds[i] = feat_inds
  PATHWAY_TO_SUPPORT = lapl_to_feat_inds
  LAPLACIANS = Ls
  for pathway_id, L in enumerate(LAPLACIANS):
    support = PATHWAY_TO_SUPPORT[pathway_id]
    L_normal = normalize_laplacian(L, support)
    NORMALIZED_LAPLACIANS.append(L_normal)

  # track which Laplacian/pathway are candidates for each latent vector
  # initially, all pathways are candidates
  latent_to_pathway_data = init_latent_to_pathway_data(k_latent, Ls)

  converged = False
  candidates_remain = True
  obj = math.inf
  obj_data = {}
  prev_obj = math.inf
  i = 0
  k_to_lapl_ind = {}
  # we make jumps in the search space when considering different pathways so some local optima 
  # may not be the best we have seen so far, track the best we have seen
  best_dict = {} 
  best_dict['obj_data'] = {}
  best_dict['obj_data']['obj'] = np.Inf
  do_multinomial = True # TODO expose as parameter
  while (i < max_iter) and (candidates_remain or not converged):
    # update active Laplacian for each latent vector every <modulus> iterations
    for k in range(k_latent):
      lapl_inds = list(map(lambda x: x[0], latent_to_pathway_data[k]))
      if do_multinomial:
        # sample current pathway where high scoring pathways are more likely to be selected
        lapl_scores = np.array(list(map(lambda x: x[1], latent_to_pathway_data[k])))
        lapl_prob = lapl_scores / np.sum(lapl_scores)
        multinomial_sample = multinomial.rvs(1, lapl_prob)
        sample_ind = np.where(multinomial_sample != 0)[0][0]
        lapl_ind = lapl_inds[sample_ind]
      else:
        # sample current pathway uniformly at random from candidate pathways
        lapl_ind = random.choice(lapl_inds)
      k_to_lapl_ind[k] = lapl_ind
    k_to_W, k_to_D, k_to_L, k_to_feat_inds = map_k_to_lapls(k_to_lapl_ind, Ws, Ds, Ls, lapl_to_feat_inds)

    if verbose:
      print('--------------------------------------------')
      print('Latent/Pathway association at this iteration')
      print('--------------------------------------------')
      for k, lapl_ind in k_to_lapl_ind.items():
        print(k, lapl_ind)
      print('--------------------------------------------')
    if tradeoff is None:
      U, V, obj_data = nmf_manifold_vec_update_normal(X, U, V, k_to_W, k_to_D, k_to_L, k_to_feat_inds, n_steps=modulus, i=i, norm_X=norm_X, gamma=gamma, delta=delta, verbose=verbose)
    else:
      U, V, obj_data, gamma, delta = nmf_manifold_vec_update_tradeoff(X, U, V, k_to_W, k_to_D, k_to_L, k_to_feat_inds, n_steps=modulus, i=i, norm_X=norm_X, tradeoff=tradeoff, gamma=gamma, delta=delta, verbose=verbose)
    i += modulus
    
    # track best
    # TODO need to track latent factor mapping too
    if obj_data['obj'] < best_dict['obj_data']['obj']:
      best_dict['U'] = U
      best_dict['V'] = V
      best_dict['obj_data'] = obj_data

    # after <modulus> updates, restrict candidates
    # dont restrict any further if the number of distinct Laplacians remaining is equal to <k_latent>
    if candidates_remain:
      if (count_distinct_pathways(latent_to_pathway_data) <= k_latent):
        # if this condition is met, force each latent factor to target different Laplacians
        latent_to_pathway_data = force_distinct_lapls(V, Ls, latent_to_pathway_data, k_to_feat_inds, gamma, delta)
        if verbose:
          print_latent_to_pathway_data(latent_to_pathway_data)
        candidates_remain = False
      else:
        sys.stderr.write('Before restrict: ' + str(datetime.datetime.now()) + '\n')
        latent_to_pathway_data = restrict(V, Ls, latent_to_pathway_data, lapl_to_feat_inds)
        sys.stderr.write('After restrict: ' + str(datetime.datetime.now()) + '\n')
        candidates_remain = False
        for k,v in latent_to_pathway_data.items():
          if(len(v) > 1):
            candidates_remain = True
        if verbose:
          print_latent_to_pathway_data(latent_to_pathway_data)

    prev_obj = obj
    obj = obj_data['obj']
    converged = (abs(obj - prev_obj)) / obj < tol
  # end while

  # replace converged result with the best observed
  if best_dict['obj_data']['obj'] < obj_data['obj']:
    print('Local optima at convergence (or after max iterations) is not the best among all iterates; returning best instead')
    U = best_dict['U']
    V = best_dict['V']
    obj_data = best_dict['obj_data']

  # TODO this preserves reconstruction error part of obj but not the other parts, report obj after rescale
  # rescale 
  # alpha X = U * V^T
  # X = U * (1/alpha * V^T)
  V = (1 / alpha) * V

  obj_data['latent_to_pathway_data'] = latent_to_pathway_data

  return U, V, obj_data

def check_header(fpath, delim):
  # check if there is a header in the data file
  has_header = False
  with open(fpath, 'r') as fh:
    for line in fh:
      line = line.rstrip()
      words = line.split(delim)
      for word in words:
        try:
          float(word)
        except ValueError as err:
          has_header=True
          break
      break
  return has_header

def check_row_names(fpath, delim, has_header):
  # check if there are row names in the data file
  has_row_names = False
  with open(fpath, 'r') as fh:
    data_line = None
    i = 0
    target_line = 1
    if has_header:
      target_line = 2
      
    for line in fh:
      i += 1
      data_line = line.rstrip()
      if i >= target_line:
        break

    words = data_line.split(delim)
    try:
      float(words[0])
    except ValueError as err:
      has_row_names = True
  return has_row_names

def parse_pathways(manifold_fps, node_attribute="name"):
  """
  Parse graphml files used for pathways. Restructure networkx objects so that gene names are 
  used as node identifiers. By default, parse the graphml file at the end of this docstring 
  so that "ENSG11486" is a node identifier rather than "11485".

    <graphml ...>
      <key attr.name="name" attr.type="string" for="node" id="d1" />
      <graph edgedefault="undirected">
        <node id="11485">
          <data key="d1">ENSG11486</data>
        </node>
        ...
      </graph>
    </graphml>

  Parameters
  ----------
  manifold_fps : list of str
    list of filepaths to pathway files (in the graphml file format)

  Returns
  -------
  G_fp_pairs : list of networkx.Graph
    parsed pathways

  TODO
  ----
  implement behavior so that "11485" could be used as a node identifier
  """
  G_fp_pairs = []
  for fp in manifold_fps:
    G = nx.read_graphml(fp).to_undirected()
    G = prmf.relabel_nodes(G, node_attribute)
    G_fp_pairs.append((G,fp))
  return G_fp_pairs

def main():
  parser = argparse.ArgumentParser(description="""
Python implementation of Pathway-Regularized NMF.

Solve an optimization problem of the form
  min ||X - UV^T|| + 
    gamma * sum_k min_i V[:,k]^T Ls[i] V[:,k] + 
    delta * sum_k sum_{i | i in G_k} 1 / V[i,k] + 
    ||U||_F^2

where Ls[i] is the Laplacian matrix associated with Gs[i],
G_k is the manifold associated with latent factor k
X has shape (n_obs, n_features),
U has shape (n_obs, n_latent),
V has shape (n_feature, n_latent)

References
----------
Cai 2008. Non-negative Matrix Factorization on Manifold
""", formatter_class=RawTextHelpFormatter)
  prmf.add_prmf_arguments(parser)
  args = parser.parse_args()
  OUTDIR = args.outdir

  # tradeoff, gamma, and delta
  tradeoff = args.tradeoff
  if tradeoff == -1:
    tradeoff = None

  # TODO update gamma default

  manifold_fps = []
  if args.manifolds is None and args.manifolds_file is None:
    sys.stderr.write("Exactly one of --manifolds or --manifolds-file is required.\n")
    sys.exit(22)
  elif args.manifolds is None and args.manifolds_file is not None:
    with open(args.manifolds_file, 'r') as fh:
      for line in fh:
        line = line.rstrip()
        manifold_fps.append(line)
  elif args.manifolds is not None and args.manifolds_file is None:
    manifold_fps = args.manifolds
  else:
    sys.stderr.write("Exactly one of --manifolds or --manifolds-file is required.\n")
    sys.exit(23)
  G_fp_pairs = parse_pathways(manifold_fps)
  fp_to_G = {}
  for G, fp in G_fp_pairs:
    fp_to_G[fp] = G
  Gs = list(map(lambda x: x[0], G_fp_pairs))

  # TODO warn if --node-attribute is not found

  if args.seed is not None:
    seed = int(args.seed)
    np.random.seed(seed)
    random.seed(seed)

  has_header = check_header(args.data, args.delimiter)
  has_row_names = check_row_names(args.data, args.delimiter, has_header)
  # load data
  X = None
  # pd.read_csv defaults updated by CLI arguments
  nrows = None
  if args.m_samples is not None:
    n_rows = args.m_samples
  header = 'infer'
  if not has_header:
    header = None
  index_col = None
  if has_row_names:
    index_col = 0
  X = pd.read_csv(args.data, sep=args.delimiter, header=header, nrows=nrows, index_col=index_col)
  samples = None
  if has_row_names:
    samples = list(X.index)
  
  # transpose data if desired
  m, n = X.shape
  if args.high_dimensional:
    if m > n:
      X = X.transpose()
  else:
    if m < n:
      X = X.transpose()

  # finalize data prep for nmf_pathway:
  # parse nodelist if provided or infer it from X as a dataframe
  # convert data frame to numpy
  nodelist = None
  if args.nodelist is not None:
    nodelist = prmf.parse_nodelist(open(args.nodelist))
    X = X.to_numpy()
  else:
    if has_header:
      # use the header to construct a nodelist
      nodelist = list(X.columns)
      nodelist_set = set(nodelist)
      for G in Gs:
        for node in G:
          if node not in nodelist_set:
            nodelist.append(node)
            nodelist_set.add(node)

      X = prmf.embed_arr(nodelist, list(X.columns), X.to_numpy())
    else:
      sys.stderr.write("--nodelist is not provided and there is no header in <--data>\n")
      sys.exit(25)

  # check node identifiers in G against nodelist
  # TODO rework this test for inferred nodelist
  nodelist_set = set(nodelist)
  G_index_to_frac = {}
  all_zero = True
  for i,G in enumerate(Gs):
    count = 0
    for node in G.nodes_iter():
      if node in nodelist_set:
        count += 1
    frac = count / G.order()
    G_index_to_frac[i] = frac
    if count != 0:
      all_zero = False
  if all_zero:
    sys.stderr.write("Invalid manifolds. Check that the node identifiers of the manifolds are present in the nodelist. Try setting --node-attribute if the node identifier is in a graphml attribute rather than the XML node attribute 'id'\n")
    sys.exit(24)
  sys.stdout.write("Printing manifold node representation in nodelist:\n")
  for i, G_fp_pair in enumerate(G_fp_pairs):
    sys.stdout.write("{}: {:2.1f}%\n".format(G_fp_pair[1], G_index_to_frac[i]*100))

  U_fp = os.path.join(args.outdir, "U.csv")
  V_fp = os.path.join(args.outdir, "V.csv")
  obj_fp = os.path.join(args.outdir, "obj.txt")

  # cross validation
  # TODO use other folds
  X_test = None
  if args.cross_validation is not None:
    kf = KFold(n_splits = round(1/args.cross_validation))
    for train_index, test_index in kf.split(X):
      X_train = X[train_index]
      X_test = X[test_index]

      X = X_train
      samples = [samples[i] for i in train_index]
      break

  # normalize data if desired
  # data at this stage is assumed to be observations x features
  # normalization is done for each feature value
  # e.g. the sample with the highest read count for gene X gets the value 1 in the gene X column
  if not args.no_normalize:
    X = quantile_transform(X)

  # --manifolds-init - {{
  pathway_init_fp = os.path.join(args.outdir, 'init_pathways.txt')
  U_init = None
  V_init = None
  init_fps = []
  if args.manifolds_init is not None:
    Gs_init = list(map(lambda fp: fp_to_G[fp], args.manifolds_init))
    if len(args.manifolds_init) < args.k_latent:
      # then extend Gs_init with a random sample from the pathway population
      non_init_fps = list(set(manifold_fps) - set(args.manifolds_init))
      chosen_fps = random.sample(non_init_fps, args.k_latent - len(args.manifolds_init))
      init_fps = copy.copy(args.manifolds_init)
      for chosen_fp in chosen_fps:
        Gs_init.append(fp_to_G[chosen_fp])
        init_fps.append(chosen_fp)
    elif len(args.manifolds_init) == args.k_latent:
      # no modification to Gs_init is needed
      init_fps = args.manifolds_init
    else: # len(args.manifolds_init) > args.k_latent
      # then sample from Gs_init
      inds = np.random.choice(len(Gs_init), args.k_latent)
      Gs_init_new = []
      for ind in inds:
        Gs_init_new.append(Gs_init[ind])
        init_fps.append(args.manifolds_init[ind])
      Gs_init = Gs_init_new
    vs = []
    us = []
    for G in Gs_init:
      v, pathway_ind = pathway_to_vec(X, G, nodelist)
      v_pathway_signal = v[pathway_ind]
      u, res = nmf_init_u(X, v)
      v_new, res = nmf_init_v(X, u)
      v_new[pathway_ind] = v_pathway_signal
      vs.append(v_new)
      us.append(u)
    V_init = np.concatenate(vs, axis=1)
    U_init = np.concatenate(us, axis=1)
    sys.stdout.write("Using the following manifolds for initialization:\n{}\n".format("\n".join(init_fps)))
    # also write these to their own file
    with open(pathway_init_fp, 'w') as pathway_init_fh:
      pathway_init_fh.write("\n".join(init_fps))
  # }} - --manifolds-init

  # TODO other arguments
  U, V, obj_data = nmf_pathway(X, Gs, nodelist=nodelist, gamma=args.gamma, tradeoff=tradeoff, k_latent=args.k_latent, U_init=U_init, V_init=V_init)
  U = pd.DataFrame(U, index=samples, columns=list(map(lambda x: "LV{}".format(x), range(args.k_latent))))
  V = pd.DataFrame(V, index=nodelist, columns=list(map(lambda x: "LV{}".format(x), range(args.k_latent))))
  U.to_csv(U_fp, sep=",", index=has_row_names, quoting=csv.QUOTE_NONNUMERIC)
  V.to_csv(V_fp, sep=",", index=True, quoting=csv.QUOTE_NONNUMERIC)

  # cross validation
  if args.cross_validation is not None:
    normalized_test_errors = prmf.measure_cv_performance(V, X_test)
    avg_normalized_test_error = np.mean(normalized_test_errors)
    error_fp = os.path.join(args.outdir, 'test_error.csv')
    np.savetxt(error_fp, normalized_test_errors, delimiter=",")
    obj_data['average_normalized_test_error'] = avg_normalized_test_error

  with open(obj_fp, 'w') as obj_fh:
    latent_to_pathway_data = obj_data.pop('latent_to_pathway_data', {})
    for k,v in obj_data.items():
      obj_fh.write("{} = {:0.5f}\n".format(k,v))

    # write which manifold file was used for each latent factor
    ks = sorted(latent_to_pathway_data.keys())
    for k in ks:
      lapl_inds = list(map(lambda x: x[0], latent_to_pathway_data[k]))
      # TODO pick first, assumes convergence
      lapl_ind = lapl_inds[0]
      G, fp = G_fp_pairs[lapl_ind]
      obj_fh.write("{} -> {}\n".format(k, fp))

if __name__ == "__main__":
  main()
