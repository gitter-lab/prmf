import sys
import math
import os, os.path
import re
import itertools as it
import numpy as np
import networkx as nx
import sklearn
import sklearn.metrics
import scipy.sparse as sp
from scipy.stats import gamma
from scipy.sparse.linalg import norm

class FactorLibException(Exception):
  pass

# TODO copied from ppi.irefindex
def transform_nearest_neighbors(G, k=10, attr='weight', attr_type='similarity'):
  """
  Remove edges in G. For each node, keep at most 10 of its neighbors. Kept neighbors have the 
  highest 'weight' edge attribute. Returns a copy of G and does not modify G.

  Parameters
  ----------
  G : nx.Graph

  k : int
    default 10

  attr : str
    default 'weight'

  attr_type : str
    'similarity': nearest neighbors are those that are most similar: keep neighbors with highest edge weight
    'distance': nearest neighbers are those that are least distant from each other: keep neighbors with least edge weight

  Returns
  -----
  H : nx.Graph
    k-nearest neighbor graph
    note that a vertex v may have more than k edges incident to it if its k'>k neighbor has v in its k-nearest neighbors
  """
  attr_types = ['similarity', 'distance']
  if(attr_type not in attr_types):
    raise ArgumentError("attr_type: {} should be one of {}".format(attr_type, ",".join(attr_types)))
  H = nx.Graph() # TODO assumes undirected
  for node in G.nodes_iter():
    inc_edges = G.edges([node])
    edge_weight_pairs = []
    for edge in inc_edges:
      weight = G[edge[0]][edge[1]][attr]
      pair = (edge, weight)
      edge_weight_pairs.append(pair)
    if(attr_type == 'similarity'):
      edge_weight_pairs = sorted(edge_weight_pairs, key=lambda x: x[1], reverse=True)
    else:
      edge_weight_pairs = sorted(edge_weight_pairs, key=lambda x: x[1])
    for pair in edge_weight_pairs[:k]:
      edge = pair[0]
      weight = pair[1]
      H.add_edge(edge[0], edge[1], {'weight': weight})
  return H

def rbf_similarity(stdev, x, y):
  return math.exp(-(np.abs(x - y) / stdev) ** 2)

def combine_graphs(Gs, weights):
  """
  Combine graphs into a single graph by constructing a weighted sum of their edge weights

  G_1 := (V_1, E_1) 
  G_2 := (V_2, E_2)
  where E_i is a set of tuples (u,v,w)

  G := (V_1 \cup V_2, E)
  TODO
  """
  def get_edge_weight(G, u, v):
    # assumes edge is already known to be present in G
    weight = G[u][v].get('weight')
    if weight is None:
      weight = 1
    return weight

  total_weight = 0.0
  for weight in weights:
    total_weight += weight

  G_rv = nx.Graph()
  for i, G in enumerate(Gs):
    graph_weight = weights[i]
    for edge in G.edges_iter():
      weight = get_edge_weight(G, *edge)

      if G_rv.has_edge(*edge):
        weight_rv = G_rv[edge[0]][edge[1]]['weight']
        G_rv[edge[0]][edge[1]]['weight'] = weight_rv + graph_weight * weight
      else:
        G_rv.add_edge(edge[0], edge[1], {'weight': graph_weight * weight})
    
  return G_rv

def diffusion(M, adj, alpha=0.7, tol=10e-6):  # TODO equation, M, alpha
    """Network propagation iterative process

    Iterative algorithm for apply propagation using random walk on a network:
        Initialize::
            X1 = M

        Repeat::
            X2 = alpha * X1.A + (1-alpha) * M
            X1 = X2

        Until::
            norm(X2-X1) < tol

        Where::
            A : degree-normalized adjacency matrix

    Parameters
    ----------
    M : sparse matrix
        Data matrix to be diffused.

    adj : sparse matrix
        Adjacency matrice.

    alpha : float, default: 0.7
        Diffusion/propagation factor with 0 <= alpha <= 1.
        For alpha = 0 : no diffusion.
        For alpha = 1 :

    tol : float, default: 10e-6
        Convergence threshold.

    Returns
    -------
    X2 : sparse matrix
        Smoothed matrix.

    Notes
    -----
    Copied from the stratipy Python library
    """
    n = adj.shape[0]
    adj = adj+sp.eye(n)

    d = sp.dia_matrix((np.array(adj.sum(axis=0))**-1, [0]), shape=(n,  n))
    A = adj.dot(d)

    X1 = M
    X2 = alpha * X1.dot(A) + (1-alpha) * M
    i = 0
    while norm(X2-X1) > tol:
        X1 = X2
        X2 = alpha * X1.dot(A) + (1-alpha) * M
        i += 1
    return X2

def parse_ws_delim(fh):
  return parse_nodelist(fh)

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

def prepare_nodelist(lists):
  """
  Combine lists of nodes into a single master list for use by nmf_pathway and other functions that
  accept a nodelist as input
  """
  set_all = set()
  for ll in lists:
    for item in ll:
      set_all.add(item)
  rv = sorted(set_all)
  return rv

def parse_nodelist(fh):
  """
  Return a list of node identifiers which maps the list index to the identifier
  """
  rv = []
  for line in fh:
    line = line.rstrip()
    words = line.split()
    for word in words:
      rv.append(word)
  return rv

def parse_pathways(pathways_dir):
  """
  Parse .graphml files in <pathways_dir>

  Returns
  -------
  Gs : list of nx.DiGraph
    parsed pathways in <pathways_dir>
  """
  graphmls = []
  for fname in os.listdir(pathways_dir):
    basename, ext = os.path.splitext(fname)
    if ext == ".graphml":
      graphmls.append(os.path.join(pathways_dir, fname))
  Gs = list(map(lambda x: nx.read_graphml(x), graphmls))
  return Gs

def embed_ids(all_ids, ids):
  """
  Construct a binary vector of length len(<all_ids>) with a 1 in the positions
  associated with each id in <ids>
  """
  row_vec = np.zeros(len(all_ids))
  missing = []
  for id in ids:
    try:
      index = all_ids.index(id)
      row_vec[index] = 1
    except ValueError:
      missing.append(id)

  row_vec_sparse = sp.csc_matrix(row_vec)
  return (row_vec_sparse, missing)

def embed_arr(all_col_names, some_col_names, arr):
  """
  Embed <arr> with len(some_col_names) number of columns in a larger array with len(all_col_names)
  number of columns. The larger array is constructed and returned.
  """
  m,n = arr.shape
  if len(some_col_names) != n:
    raise FactorLibException('some_col_names != #columns of arr: {} != {}'.format(some_col_names, n))

  n2 = len(all_col_names)
  rv = np.zeros((m,n2))

  all_name_to_ind = {}
  for i,name in enumerate(all_col_names):
    all_name_to_ind[name] = i

  for i in range(m):
    for j, name in enumerate(some_col_names):
      j2 = all_name_to_ind[name]
      rv[i,j2] = arr[i,j]

  return rv

def get_adj_mat(G):
  """Represent ppi network as adjacency matrix

  Parameters
  ----------
  G : networkx graph
    ppi network, see get_ppi()

  Returns
  -------
  adj : square sparse scipy matrix
    (i,j) has a 1 if there is an interaction reported by irefindex

  ids : list
    same length as adj, ith index contains irefindex unique identifier for gene whose interactions are
    reported in the ith row of adj
  """
  ids = G.nodes()
  adj = nx.to_scipy_sparse_matrix(G, nodelist=ids, dtype=bool)

  return adj, ids

def weighted_union(Gs, weight_attr='weight'):
  """
  Construct a weighted graph that is the union of graphs <Gs>.
  For each edge in a network in <Gs>, use <weight_attr> as the edge weight or 1 if that attribute is not present.
  Each edge in G_rv has its <weight_attr> set to the average weight of <Gs> that the edge appears in.
  That is, networks without the edge do not pull down the average toward 0.

  Parameters
  ----------
  Gs : iterable of nx.Graph

  weight_attr : str
    name of edge attribute to use as weight

  Returns
  -------
  G_rv : nx.Graph
  """
  G_rv = nx.Graph()
  edge_to_count = {}

  def get_edge_weight(G, edge):
    weight = G[edge[0]][edge[1]].get(weight_attr)
    if weight is None:
      weight = 1
    return weight

  for G in Gs:
    for edge in G.edges_iter():
      weight = get_edge_weight(G, edge)
      prev_weight = 0
      if G_rv.has_edge(*edge):
        prev_weight = get_edge_weight(G_rv, edge)
      weight_sum = prev_weight + weight
      G_rv.add_edge(edge[0], edge[1], {weight_attr: weight_sum})
      if edge in edge_to_count:
        edge_to_count[edge] += 1
      else:
        edge_to_count[edge] = 1

    for node in G.nodes_iter():
      if node not in G_rv:
        G_rv.add_node(node)

  for edge, count in edge_to_count.items():
    weight = get_edge_weight(G_rv, edge)
    avg_weight = float(weight) / count
    G_rv[edge[0]][edge[1]][weight_attr] = avg_weight

  return G_rv

def mat_to_bipartite(mat):
  G = nx.Graph()
  row_node_ids = []
  col_node_ids = []
  n_rows, n_cols = mat.shape
  for i in range(n_rows):
    node_i = "r{}".format(i)
    row_node_ids.append(node_i)
    G.add_node(node_i, bipartite=0)

  for j in range(n_cols):
    node_j = "c{}".format(j)
    col_node_ids.append(node_j)
    G.add_node(node_j, bipartite=1)

  for i in range(n_rows):
    for j in range(n_cols):
      G.add_edge("r{}".format(i), "c{}".format(j), {'weight': mat[i,j]})
  return G, row_node_ids, col_node_ids

def to_unit_vector(vec):
  rv = None
  norm = np.linalg.norm(vec)
  if(norm == 0):
    rv = vec
  else:
    rv = vec / norm
  return rv

def nodelists_to_mat(sub_lists, uni_list):
  """Represent gene nodelist as a matrix with shape (len(uni_list), len(sub_lists))

  Parameters
  ----------
  sub_lists : list of str

  uni_list : list of str

  Returns
  -------
  mat : np.array
    matrix where each column is associated with a nodelist in <sub_lists> and has a 1 in each
    row where that gene is in the node list

  TODO
  ----
  For a single pathway, this representation has n_genes_in_pathway "units of signal" while
  the sampled signal has units <= n_genes_in_pathway (because a subset is sampled)
  """
  node_to_index = {}
  for i, node in enumerate(uni_list):
    node_to_index[node] = i

  n_pathways = len(sub_lists)
  mat = np.zeros((len(uni_list), n_pathways))
  for pathway_ind in range(len(sub_lists)):
    pathway = sub_lists[pathway_ind]
    for node in pathway:
      index = node_to_index[node]
      if index is None:
        sys.stderr.write("Unrecognized gene ID: {}\n".format(node))
      else:
        mat[index,pathway_ind] = 1
  return mat

def normalize_num_pathways(W_mat, pathways_mat):
  """
  Prior to matching latent factors to pathways, this function MAY be invoked to transform
  W_mat or pathways_mat (whichever is smaller) to make it so the number of latent factors
  is the same as the number of pathways by adding a binary vector of all zeros (the empty pathway).

  Parameters
  ----------
  W_mat : np.array
    array with latent factors as columns and genes as rows with gene loadings into latent factors in cells

  pathways_mat : np.array
    array with pathways as columns and genes as rows with a 1 or 0 in each cell indicating whether that
    gene is present in the pathway or not

  Returns
  -------
  W_mat : np.array
    See description

  pathways_mat : np.array
    See description
  """
  # reshape 1D vectors if needed
  shp = W_mat.shape
  if(len(shp) == 1):
    W_mat = W_mat.reshape((shp[0], 1))
  shp = pathways_mat.shape
  if(len(shp) == 1):
    pathways_mat = pathways_mat.reshape((shp[0], 1))

  n_genes, n_latent_factors = W_mat.shape
  shp = pathways_mat.shape
  n_genes2, n_pathways = pathways_mat.shape
  if(n_genes != n_genes2):
    raise ValueError("{} != {}: number of genes identified by W_mat is not the same as those identified by pathways_mat".format(n_genes, n_genes2))

  arrs = [W_mat, pathways_mat]
  small_arr = -1
  large_arr = -1
  diff = 0
  if(n_latent_factors < n_pathways):
    small_arr = 0
    large_arr = 1
    diff = n_pathways - n_latent_factors
  elif(n_pathways < n_latent_factors):
    small_arr = 1
    large_arr = 0
    diff = n_latent_factors - n_pathways
  # else do nothing

  cols = np.zeros((n_genes, diff))
  arrs[small_arr] = np.concatenate((arrs[small_arr], cols), axis=1)

  return tuple(arrs)

def match(W_mat, pathways_mat):
  return match_pr(W_mat, pathways_mat)

def match_pr(W_mat, pathways_mat):
  n_genes, n_factors = W_mat.shape
  n_genes2, n_pathways = pathways_mat.shape
  if(n_genes != n_genes2):
    raise ValueError("n_genes != n_genes2 : {} != {}".format(n_genes, n_genes2))

  # prepare max weight matching by AUC
  aucs = np.zeros((n_factors, n_pathways))
  for i in range(n_factors):
    factor_vec = W_mat[:,i]
    for j in range(n_pathways):
      pathway_vec = pathways_mat[:,j]

      y_score = factor_vec
      y_true = pathway_vec
      auc = sklearn.metrics.average_precision_score(y_true, y_score)
      aucs[i,j] = auc

  # use networkx for max weight matching
  G, factor_node_ids, pathway_node_ids = mat_to_bipartite(aucs)
  matching = nx.max_weight_matching(G)

  # I don't like the return value from networkx, reorganize data
  rv = transform_matching(G, matching, factor_node_ids)
  return rv

def match_dist(W_mat, pathways_mat):
  """
  Match latent factors to pathways

  Parameters
  ----------
  W_mat : np.array

  pathways_mat : np.array

  Returns
  -------
  rv : list of tpl
    each tpl is of the form (<latent_factor_id>, <pathway_id>, <distance>)
    where W_mat latent factors are identified with names "r0", "r1", etc. (for "row" of the 
    distance matrix constructed in this function) associated with each latent factor in W_mat; 
    and pathways are identified with names "c0", "c1", etc. (for "column" of the distance 
    matrix)

  TODO
  ----
  other type of transformation of distance to similarity? RBF?
  implement version that uses hypergeometic p-values as distance measure?
  """
  n_genes, n_factors = W_mat.shape
  n_genes2, n_pathways = pathways_mat.shape
  if(n_genes != n_genes2):
    raise ValueError("n_genes != n_genes2 : {} != {}".format(n_genes, n_genes2))
  #if(n_factors != n_pathways):
  #  raise ValueError("n_factors != n_pathways: {} != {}".format(n_factors, n_pathways))

  dists = np.zeros((n_factors, n_pathways))
  max_unit_vec_dist = math.sqrt(2)
  for i in range(n_factors):
    factor_vec = W_mat[:,i]
    # scale distance by size of latent factor and pathway
    factor_vec_unit = to_unit_vector(factor_vec)
    for j in range(n_pathways):
      pathway_vec = pathways_mat[:,j]
      pathway_vec_unit = to_unit_vector(pathway_vec)
      dists[i,j] = np.linalg.norm(factor_vec_unit - pathway_vec_unit)
  dists_inv = 1 - (dists / max_unit_vec_dist)

  # use networkx for max weight matching
  G, factor_node_ids, pathway_node_ids = mat_to_bipartite(dists_inv)
  matching = nx.max_weight_matching(G)

  # I don't like the return value from networkx, reorganize data
  rv = transform_matching(G, matching, factor_node_ids)
  return rv

def transform_matching(G, matching, factor_node_ids):
  rv = []
  for factor_id in factor_node_ids:
    # not all nodes used in matching necessarily
    if factor_id in matching:
      pathway_id = matching[factor_id]
      weight = G[factor_id][pathway_id]['weight']
      tpl = (factor_id, pathway_id, weight)
      rv.append(tpl)
    else:
      # TODO temporary: why are they not included?
      # usually not included due to different numbers of latent factors and pathways
      weights = []
      for neighbor in G.neighbors(factor_id):
        weights.append(G[factor_id][neighbor]['weight'])
      sys.stderr.write("Factor {} not included in matching; weights: ".format(factor_id) + " ; ".join(map(str, weights)) + "\n")
  return rv

def sample(G, n, t=3, alpha_portion=0.9, seed=None, node_order=None):
  """
  Sample <n> random values according to a distribution specified by G.
  Each value is a len(node_order) size vector of positive reals.

  If a node in G has no parents, its value is drawn according to Gamma(t^alpha_portion, t^(1-alpha_portion)).
  t is a reparameterization so that t is a location parameter.
  The default alpha_portion is chosen to be high so that the distribution is less uniform and more peaked around that location.

  Nodes which are present in <node_order> but not <G> are also have values drawn according to the aforementioned distribution.

  A node Y which has parents X_1, ..., X_k is specified by Y = \sum_k Gamma(X_k^alpha_portion, X_k^(1-alpha_portion))
  Y also has a Gamma distribution but it is not easily expressed. For Y = \sum_k Gamma(alpha_k, beta),
  Y ~ Gamma(\sum_k alpha_k, beta). The Y here is different but close.
  """
  loc = 0
  # TODO assumes connected
  if not G.is_directed():
    raise FactorLibException('G must be directed')
  if seed is not None:
    np.random.seed(seed)
  if node_order is None:
    node_order = G.nodes()

  m = len(node_order)
  rv = np.zeros((m,n))
  node_to_ind = {}
  for i, node in enumerate(node_order):
    node_to_ind[node] = i

  roots = [n for n,d in G.in_degree().items() if d==0]
  for j in range(n):
    for root in roots:
      root_ind = node_to_ind[root]
      alpha = math.pow(t, alpha_portion)
      beta = math.pow(t, 1 - alpha_portion)
      rv[root_ind,j] = gamma.rvs(alpha, loc, beta)
      for source, target in nx.dfs_edges(G, root):
        source_ind = node_to_ind[source]
        target_ind = node_to_ind[target]
        alpha = math.pow(rv[source_ind,j], alpha_portion)
        beta = math.pow(rv[source_ind,j], 1 - alpha_portion)
        rv[target_ind,j] += gamma.rvs(alpha, loc, beta)

  other_nodes = set(node_order) - set(G.nodes())
  for j in range(n):
    for node in other_nodes:
      node_ind = node_to_ind[node]
      alpha = math.pow(t, alpha_portion)
      beta = math.pow(t, 1 - alpha_portion)
      rv[node_ind, j] = gamma.rvs(alpha, loc, beta)

  return rv

def compute_man_reg(w, L):
  return (L.dot(w)).dot(w)

def parse_achilles(fp):
  """
  Parse Achilles data from http://portals.broadinstitute.org/achilles/datasets/18/download

  Parameters
  ----------
  fp : str
    filepath to 'gene_dependency' file mentioned in Achilles README
  """
  data = None
  row_to_cell_line = []
  gene_names = None

  with open(fp) as fh:
    for line in fh:
      line = line.rstrip()
      gene_names = line.split(',')[1:]
      break

    for line in fh:
      ind = line.index(',')
      cell_line = line[0:ind]
      row_to_cell_line.append(cell_line)

  data = np.genfromtxt(fp, delimiter=",", skip_header=1)
  data = data[:, 1:]

  gene_names = achilles_to_hugo(gene_names)

  return data, row_to_cell_line, gene_names

def achilles_to_hugo(gene_names):
  """
  Transform gene identifiers from "HUGO (Entrez)" to "HUGO"
  """
  regexp = re.compile(r'([\w-]+)\s\(\w+\)')
  rv = []
  for i, gene_name in enumerate(gene_names):
    match_data = regexp.match(gene_name)
    if match_data is not None:
      rv.append(match_data.group(1))
    else:
      raise FactorLibException("Gene name {} at position {} does not match the specified format".format(gene_name, i))
  return rv

def filter_vec_by_graph(G, vec, nodelist):
  """
  Return a subset of vec and nodelist corresponding to the nodes defined in G
  """
  inds = []
  nodelist_sub = []
  for i, node in enumerate(nodelist):
    if node in G:
      inds.append(i)
      nodelist_sub.append(node)
  vec_sub = vec[inds]
  return vec_sub, nodelist_sub
