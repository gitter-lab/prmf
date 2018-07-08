#!/usr/bin/env python
import argparse, sys
import numpy as np
import scipy.io
import scipy.sparse as sp
import networkx as nx 
import factorlib as fl

def main():
  parser = argparse.ArgumentParser(description="""
Diffuse node scores over a network. The diffused matrix is of shape (n_nodes, n_gene_lists) === (n_feature x n_obs).
""")
  parser.add_argument("--network", type=str, help="graphml network file to run diffusion on", required=True)
  parser.add_argument("--nodelist", type=argparse.FileType("r"), required=True,
    help="Association between gene identifier and matrix index provided as a whitespace delimited list")
  parser.add_argument("--gene-lists", nargs="+", help="one or more files with an node identifiers on each line")
  parser.add_argument("--gene-csv", help="One csv file with genes along columns and observations along rows; must contain column names but not row names")
  parser.add_argument("--diffused", "-d", type=argparse.FileType("wb"), required=True,
    help="Diffused matrix")
  parser.add_argument("--alpha", "-a", type=float, default=0.7,
    help="Diffusion rate parameter")
  parser.add_argument("--tolerance", "-t", type=float, default=10e-6,
    help="Tolerance threshold for diffusion; stop when change in diffused matrix crosses below threshold")
  parser.add_argument("--string-edge-type", default="combined_score", help="\"experimental\" for edges supported by experimental evidence only; \"combined_score\" for the entire stringdb network; default=\"combined_score\"")
  parser.add_argument("--diffused-format", type=str, default='csv', help="Either \"ampl\" or \"csv\"; default=\"csv\" which is short for MatrixMarket, a sparse matrix file format")
  args = parser.parse_args()

  # TODO ampl only right now
  # fail fast on --diffused-format
  #if args.diffused_format not in ['ampl', 'mm']:
  #  sys.stderr.write("invalid --diffused-format={}\n".format(args.diffused_format))
  #  sys.exit(22)

  if args.gene_lists is None and args.gene_csv is None:
    sys.stderr.write("Exactly one of --gene-lists or --gene-csv is required")
    sys.exit(23)

  # TODO edge confidence threshold, edge_type in other script
  G_ppi = nx.read_graphml(args.network)
  nodelist = fl.parse_nodelist(args.nodelist)

  # NOTE if G_ppi has 'weight' attribute on edges, its value is used; otherwise a value of 
  # 1 is populated in the ij entry for an edge (i, j)
  adj = nx.to_scipy_sparse_matrix(G_ppi, nodelist=nodelist, dtype=bool)

  mat = None
  if args.gene_lists is not None:
    # parse gene lists
    gene_lists = []
    for gene_path in args.gene_lists:
      with open(gene_path) as fh:
        gene_lists.append(fl.parse_ws_delim(fh))

    # verify gene lists present in ppi_db
    def get_row_vec_for_gene_list(gene_list):
      row_vec, missing = fl.embed_ids(nodelist, gene_list)
      sys.stderr.write("missing {}/{} node identifiers: {}\n".format(len(missing), len(gene_list), ", ".join(missing)))
      return row_vec
    row_vecs = map(get_row_vec_for_gene_list, gene_lists)

    mat = sp.vstack(row_vecs)
  else:
    mat = sp.csc_matrix(np.genfromtxt(args.gene_csv, delimiter=","))

  # do diffusion
  smoothed_mat = fl.diffusion(mat, adj, alpha=args.alpha, tol=args.tolerance)

  # write results
  if args.diffused_format == "ampl":
    # TODO does this work with 'wb'?
    fl.ampl_write_sparse_arr(smoothed_mat, args.diffused, len(nodelist))
  else:
    np.savetxt(args.diffused, smoothed_mat.todense(), delimiter=",")

if __name__ == "__main__":
  main()
