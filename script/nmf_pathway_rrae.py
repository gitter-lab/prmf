#!/usr/bin/env python
import argparse, sys
import os, os.path
import networkx as nx
from prmf import script_utils
from prmf import prmf_args

def main():
  parser = argparse.ArgumentParser(description="""
Run nmf_pathway.py with different _r_andom _r_estarts _a_nd _e_dges.

Pass through arguments to nmf_pathway.py except for --condor, --manifolds-init, --n-runs, and --outdir

TODO
----
Capture --manifolds and --manifolds-file and randomize the edges before passing to nmf_pathway.py
""")
  parser.add_argument("--n-runs", type=int, default=2, help="Number of random restarts")
  parser.add_argument("--condor", action='store_true', help='If flag is provided, submit jobs to Condor rather than running them on this machine')
  prmf_args.add_prmf_arguments(parser)
  args = parser.parse_args()

  job_graph = nx.DiGraph()
  job_id = 0
  args_dict = vars(args)
  outdir = args_dict.pop('outdir')
  condor = args_dict.pop('condor')


  # TODO what if manifolds-file is provided instead?
  manifolds = args_dict.pop('manifolds')
  args_dict.pop('manifolds_file')

  args_dict.pop('manifolds_init')
  n_runs = args_dict.pop('n_runs')
  for i in range(n_runs):
    run_outdir = os.path.join(outdir, "run{}".format(i))
    os.mkdir(run_outdir)

    random_pathway_dir = os.path.join(run_outdir, 'random_pathways')
    os.mkdir(random_pathway_dir)
    attrs = {
      'exe': 'randomize_network.R',
      'args': ['--infiles'] + manifolds + ['--outdir', random_pathway_dir],
      'out': os.path.join(random_pathway_dir, 'randomize_network.out'),
      'err': os.path.join(random_pathway_dir, 'randomize_network.err')
    }
    randomize_network_job_id = job_id
    job_graph.add_node(randomize_network_job_id, attrs)
    job_id += 1
    pathways_file = os.path.join(random_pathway_dir, 'pathways_file.txt')

    args_list = script_utils.args_to_list(args_dict)
    args_list = args_list + ['--outdir', run_outdir, '--manifolds-init', '--manifolds-file', pathways_file]
    attrs = {
      'exe': 'nmf_pathway.py',
      'args': args_list,
      'out': os.path.join(run_outdir, 'nmf_pathway.out'),
      'err': os.path.join(run_outdir, 'nmf_pathway.err')
    }
    nmf_pathway_job_id = job_id
    job_graph.add_node(nmf_pathway_job_id, attrs)
    job_graph.add_edge(randomize_network_job_id, nmf_pathway_job_id)
    job_id += 1
  job_ids = script_utils.run_digraph(outdir, job_graph, condor=condor)

if __name__ == "__main__":
  main()
