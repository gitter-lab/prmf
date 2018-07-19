#!/usr/bin/env python
import argparse, sys
import numpy as np
import scipy.sparse as sp
import networkx as nx
import factorlib as fl
import math
import os, os.path
import random
import psutil
import datetime
import pdb
np.seterr(divide='raise')
EPSILON = np.finfo(np.float32).eps

# TODO remove
OUTDIR = None

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

def restrict(V, Ls, ind_to_lapls, lapl_to_feat_inds):
  """
  For each column vector index pair (v,i) : v = V[:,i], find a subset of the Laplacians 
  Ls := ind_to_lapls[i] which have the property that the nodes in the graph associated
  with L_j are relatively overrepresented in the mass of v:
    
    TODO (relative change)
    \sum_j w[j : L_{jj} != 1] / sum_j w_j > #{j : L_{jj} != 1} / #{j}

  Parameters
  ----------
  V : np.array

  ind_to_lapls : dict<int, int>
    key is a column index for V
    value is a Laplacian index in Ls

  Returns
  -------
  rv : dict<int, int>
  """
  n_feature, k_latent = V.shape
  rv = {}
  for k in range(k_latent):
    L_inds = ind_to_lapls[k]
    if(len(L_inds) > 1):
      v = V[:,k] 
      scores = np.zeros(len(L_inds))
      for i, L_ind in enumerate(L_inds):
        nz_inds = lapl_to_feat_inds[L_ind]
        scores[i] = (np.sum(v[nz_inds]) / len(nz_inds)) / np.sum(v)
      score_inds = np.where(scores > np.percentile(scores,19.9))[0]

      lapl_inds = []
      for score_ind in score_inds:
        lapl_ind = L_inds[score_ind]
        lapl_inds.append(lapl_ind)
      rv[k] = lapl_inds
    else:
      # otherwise converged to final Laplacian
      rv[k] = L_inds
  return rv

def init_ind_to_lapls(k_latent, Ls):
  rv = {}
  for i in range(k_latent):
    rv[i] = list(range(len(Ls)))
  return rv

# TODO rename
def count_distinct_lapls(ind_to_lapls):
  rv = np.inf
  for ind, lapls in ind_to_lapls.items():
    rv = min(rv, len(lapls))
  return rv

# TODO rename ind_to_lapls to make clear that it is mapping to Laplacian indexes and not the Laplacian matrices
def force_distinct_lapls(V, Ls, ind_to_lapls):
  G = nx.Graph()
  for k, lapls in ind_to_lapls.items():
    for lapl in lapls:
      L = Ls[lapl]
      denom = L.dot(V[:,k]).dot(V[:,k])
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
    ind_to_lapls[k_node] = [l_node]
  return ind_to_lapls

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

def nmf_manifold_vec_obj(X, U, V, k_to_L, k_to_feat_inds, gamma=1, delta=1):
  obj_recon = np.linalg.norm(X - U.dot(V.transpose()))

  obj_manifold = 0.0
  for k, L in k_to_L.items():
    obj_manifold += L.dot(V[:,k]).dot(V[:,k])
  obj_manifold = obj_manifold

  obj_ign = 0.0
  for k, nz_inds in k_to_feat_inds.items():
    obj_ign = np.sum(np.power(V[nz_inds,k] + 1, -1))

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

def nmf_manifold_vec_update(X, U, V, k_to_W, k_to_D, k_to_L, k_to_feat_inds, n_steps=10, gamma=1.0, delta=10.0, i=0, verbose=True, norm_X=None, tradeoff=0.5):
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

    obj_data = nmf_manifold_vec_obj(X, U, V, k_to_L, k_to_feat_inds, gamma=gamma)
    if(verbose):
      print(i+n_step+1, obj_data['obj'])
      print(obj_data)

  return U, V, obj_data

def nmf_manifold_vec_update_tradeoff_cur(X, U, V, k_to_W, k_to_D, k_to_L, k_to_feat_inds, n_steps=10, i=0, verbose=True, norm_X=None, tradeoff=0.5):
  """
  See nmf_manifold_vec_update; this version tries to update the _current_ gradient descent step size to match the tradeoff parameter

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
      V_up_num_man[:,k] = W.dot(V[:,k])
      V_up_denom_man[:,k] = D.dot(V[:,k])

      nz_inds = k_to_feat_inds[k]
      V_up_num_ign[nz_inds,k] = np.power(V[nz_inds,k] + 1, -2)

    # transform each gradient direction to a unit vector
    # TODO math here doesnt separately nicely, these expressions ignore the correctness to make a separation into recon and man parts
    V_up_recon = np.divide(V_up_num_recon, V_up_denom_recon, out=np.ones_like(V_up_num_recon), where=V_up_denom_recon!=0)
    norm_recon = np.linalg.norm(V_up_recon)
    V_up_recon_unit = V_up_recon / norm_recon
    V_up_num_man_ign = V_up_num_man + V_up_num_ign
    V_up_man = np.divide(V_up_num_man_ign, V_up_denom_man, out=np.ones_like(V_up_num_man_ign), where=V_up_denom_man!=0)
    norm_man = np.linalg.norm(V_up_man)
    V_up_man_unit = V_up_man / norm_man
    V_up_dir = tradeoff * V_up_recon_unit + (1 - tradeoff) * V_up_man_unit
    V_up_mag = np.sqrt(norm_man + norm_recon)
    #V_up_mag = norm_recon + norm_man
    V_up = V_up_mag * V_up_dir
    V = np.multiply(V, V_up)
    V[V < EPSILON] = EPSILON

    obj_data = nmf_manifold_vec_obj(X, U, V, k_to_L, k_to_feat_inds)
    if(verbose):
      print(i+n_step+1, obj_data['obj'])
      print(obj_data)

  return U, V, obj_data

def nmf_manifold_vec_update_tradeoff(X, U, V, k_to_W, k_to_D, k_to_L, k_to_feat_inds, n_steps=10, i=0, verbose=True, norm_X=None, tradeoff=0.5, gamma=1, delta=1):
  """
  See nmf_manifold_vec_update; this version tries to update the _next_ gradient descent step size to match the tradeoff parameter

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
    # update gamma and delta such that at the next iteration the reconstruction error and the
    # manifold regularization penalty contribute approximately equal values to the objective function
    gamma = obj_data['recon'] / obj_data['manifold']
    delta = gamma

    if(verbose):
      print(i+n_step+1, obj_data['obj'])
      print(obj_data)

  return U, V, obj_data, gamma, delta

def nmf_pathway(X, Gs, gamma=1.0, delta=1.0, tradeoff=None, k_latent=6, tol=1e-4, max_iter=1000, nodelist=None, modulus=10):
  """
  Solve an optimization problem of the form
    min ||X - UV^T|| + gamma * sum_k min_i V[:,k]^T Ls[i] V[:,k] + ||U||_F^2

  where Ls[i] is the Laplacian matrix associated with Gs[i]
  
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
    number of latent factors to decompose X into: for X \approx UV^T, U is shape (m_obs, k_latent)
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
  # rescale to prevent underflow errors
  #alpha = 1 / np.min(X[X != 0])
  alpha = 1
  X = alpha * X
  norm_X = np.linalg.norm(X)
  print('norm(X) = {}'.format(norm_X))
  m,n = X.shape
  # 1 - [0,1) \in (0,1] ; need strictly positive
  # TODO other initialization strategies
  U = 3 * (1 - np.random.rand(m,k_latent))
  V = 3 * (1 - np.random.rand(n,k_latent))

  node_to_index = {}
  for i, node in enumerate(nodelist):
    node_to_index[node] = i

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

  ind_to_lapls = init_ind_to_lapls(k_latent, Ls)

  converged = False
  candidates_remain = True
  obj = math.inf
  obj_data = {}
  prev_obj = math.inf
  i = 0
  k_to_lapl_ind = {}
  gamma = 1.0
  delta = 1.0
  while (i < max_iter) and (candidates_remain or not converged):
    # update active Laplacian for each latent vector every <modulus> iterations
    # TODO handle errors either here or restrict step: if ind_to_lapls[k] is empty random.choice will raise IndexError
    # can become empty if all pathways score the same, for example
    for k in range(k_latent):
      lapl_ind = random.choice(ind_to_lapls[k])
      k_to_lapl_ind[k] = lapl_ind
    k_to_W, k_to_D, k_to_L, k_to_feat_inds = map_k_to_lapls(k_to_lapl_ind, Ws, Ds, Ls, lapl_to_feat_inds)

    # DEBUG print choice of lapls
    print('----')
    for k, lapl_ind in k_to_lapl_ind.items():
      print(k, lapl_ind)
    print('----')
    if tradeoff is None:
      U, V, obj_data = nmf_manifold_vec_update(X, U, V, k_to_W, k_to_D, k_to_L, k_to_feat_inds, n_steps=modulus, i=i, norm_X=norm_X, gamma=gamma, delta=delta)
    else:
      U, V, obj_data, gamma, delta = nmf_manifold_vec_update_tradeoff(X, U, V, k_to_W, k_to_D, k_to_L, k_to_feat_inds, n_steps=modulus, i=i, norm_X=norm_X, tradeoff=tradeoff, gamma=gamma, delta=delta)
    i += modulus

    # after <modulus> updates, restrict candidates
    # dont restrict any further if the number of distinct Laplacians remaining is equal to <k_latent>
    if(count_distinct_lapls(ind_to_lapls) <= k_latent):
      # if this condition is met, force each latent factor to target different Laplacians
      ind_to_lapls = force_distinct_lapls(V, Ls, ind_to_lapls)
      for k,v in ind_to_lapls.items():
        print(k,v)
      candidates_remain = False
    else:
      ind_to_lapls = restrict(V, Ls, ind_to_lapls, lapl_to_feat_inds)
      candidates_remain = False
      for k,v in ind_to_lapls.items():
        if(len(v) > 1):
          candidates_remain = True
        print(k, v)

    prev_obj = obj
    obj = obj_data['obj']
    converged = (abs(obj - prev_obj)) / obj < tol
  # end while

  # rescale 
  # alpha X = U * V^T
  # X = U * (1/alpha * V^T)
  V = (1 / alpha) * V

  obj_data['ind_to_lapls'] = ind_to_lapls

  return U, V, obj_data

def main():
  parser = argparse.ArgumentParser(description="""
Python implementation of Pathway-Regularized NMF.

Solve min ||X - UV^T|| + gamma * sum_k min_i V[:,k]^T Ls[i] V[:,k] + ||U||_F^2

where X has shape (n_obs, n_features)
  U has shape (n_obs, n_latent)
  V has shape (n_feature, n_latent)

References
----------
Cai 2008. Non-negative Matrix Factorization on Manifold
""")
  parser.add_argument("--data", type=argparse.FileType('rb'), required=True, help="n_obs x n_features matrix")
  parser.add_argument("--manifolds", nargs='+', help="graphml files to use as manifold")
  parser.add_argument("--manifolds-file", help="A file containing newline-delimited filepaths which are used as graphml files as in <manifolds>")
  parser.add_argument("--outdir", type=str, required=True, help="Directory containing results")
  parser.add_argument("--nodelist", type=argparse.FileType('r'), help="Association of node identifier to matrix indexes", required=True)
  parser.add_argument("--tolerence", type=float, default=1e-4)
  parser.add_argument("--seed", default=None)
  parser.add_argument("--gamma", default=1.0, help="Tradeoff between reconstruction error and manifold regularization term; Default = 1.0", type=float)
  parser.add_argument("--delta", default=1.0, help="Regularization parameter for penalty for ignoring manifold; Default = 1.0", type=float)
  parser.add_argument("--tradeoff", default=None, help="If set, use the previous iterations objective function components to automatically update gamma and delta for the next iteration.")
  parser.add_argument("--high-dimensional", default=True, type=bool, help="If True, ensure that <data> is of shape m x n with m < n ; otherwise ensure that X is m x n with m > n")
  args = parser.parse_args()
  OUTDIR = args.outdir

  tradeoff = None
  if args.tradeoff is not None:
    tradeoff = float(args.tradeoff)

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
  G_fp_pairs = []
  for fp in manifold_fps:
    G = nx.read_graphml(fp).to_undirected()
    G_fp_pairs.append((G,fp))
  G_fp_pairs = sorted(G_fp_pairs, key=lambda x: x[0].order())
  Gs = map(lambda x: x[0], G_fp_pairs)

  if args.seed is not None:
    seed = int(args.seed)
    np.random.seed(seed)
    random.seed(seed)

  nodelist = fl.parse_nodelist(args.nodelist)

  U_fp = os.path.join(args.outdir, "U.csv")
  V_fp = os.path.join(args.outdir, "V.csv")
  obj_fp = os.path.join(args.outdir, "obj.txt")

  X = np.genfromtxt(args.data, delimiter=",")
  m, n = X.shape
  if args.high_dimensional:
    if m > n:
      X = X.transpose()
  else:
    if m < n:
      X = X.transpose()

  # TODO other arguments
  U, V, obj_data = nmf_pathway(X, Gs, nodelist=nodelist, gamma=args.gamma, tradeoff=tradeoff)
  np.savetxt(U_fp, U, delimiter=",")
  np.savetxt(V_fp, V, delimiter=",")
  with open(obj_fp, 'w') as obj_fh:
    ind_to_lapls = obj_data.pop('ind_to_lapls', {})
    for k,v in obj_data.items():
      obj_fh.write("{} = {}\n".format(k,v))

    # write which manifold file was used for each latent factor
    ks = sorted(ind_to_lapls.keys())
    for k in ks:
      lapl_inds = ind_to_lapls[k]
      # TODO pick first, assumes convergence
      lapl_ind = lapl_inds[0]
      G, fp = G_fp_pairs[lapl_ind]
      obj_fh.write("{} -> {}\n".format(k, fp))

if __name__ == "__main__":
  main()
