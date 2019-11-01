"""
Functions used in nmf_pathway.py that are shared by other scripts
"""
import argparse

def add_prmf_arguments(parser):
  parser.add_argument("--data", type=str, required=True, help="n_obs x n_features matrix")
  parser.add_argument("--manifolds", nargs='+', help="graphml files to use as manifold. Node identifiers must appear in nodelist.")
  parser.add_argument("--manifolds-file", help="A file containing newline-delimited filepaths which are used as graphml files as in <manifolds>")
  parser.add_argument("--manifolds-init", nargs='*', help='If provided, use this list of manifolds to initialize PRMF. If given as a flag or the number of arguments is less than --k-latent, add randomly chosen manifolds to the initialization until we have --k-latent manifolds. If exactly --k-latent files are given, use those to initialize PRMF. If greater than --k-latent files are given, select a sample of size --k-latent from them.')
  parser.add_argument("--node-attribute", help="Relabel nodes in manifolds/graphs so that their node identifiers come from this node attribute. If None, node identifiers are left as is", default=None)
  parser.add_argument("--outdir", type=str, required=True, help="Directory containing results")
  parser.add_argument("--nodelist", type=str, help="Association of node identifier to matrix indexes. If not provided, will be inferred from the header in <--data> and will fail if there is no header.")
  parser.add_argument("--k-latent", "-k", default=6, help="Number of latent factors", type=int)
  parser.add_argument("--tolerence", type=float, default=1e-3)
  parser.add_argument("--seed", default=None)
  parser.add_argument("--gamma", default=1.0, help="Tradeoff between reconstruction error and manifold regularization term; Default = 1.0", type=float)
  parser.add_argument("--delta", default=1.0, help="Regularization parameter for penalty for ignoring manifold; Default = 1.0", type=float)
  parser.add_argument("--tradeoff", default=-1, type=float, help="If set, use the previous iterations objective function components to automatically update gamma and delta for the next iteration. Must be in [0,1]. Higher values place higher importance on the manifold regularization term over the reconstruction term. To disable automatic updating of gamma and delta, set tradeoff to -1. Default = -1.")
  parser.add_argument("--high-dimensional", default=True, type=bool, help="If True, ensure that <data> is of shape m x n with m < n ; otherwise ensure that X is m x n with m > n. Default = True.")
  parser.add_argument("--no-normalize", action='store_true', help="If flag is provided, don't quantile normalize the data")
  parser.add_argument("--delimiter", default=",", help="Field delimiter in <--data>")
  parser.add_argument("--m-samples", help="If provided, only use the first <--m-samples> rows in <--data>; useful for observing script behavior on a smaller dataset", type=int)
  parser.add_argument("--cross-validation", "-c", type=float, help="If provided, use the --cross-validation value as a fraction of the samples to hold out and measure model performance with")
