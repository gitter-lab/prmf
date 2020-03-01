#!/usr/bin/env python
import argparse, sys
import os, os.path
import networkx as nx
from prmf import script_utils

def main():
  parser = argparse.ArgumentParser(description=
"""
Evalute nmf_pathway.py using true pathways against randomized pathways for a "real" dataset instead of a simulated one
""")
  parser.add_argument("--rng-seed", help="Seed for random number generators", default=None)
  parser.add_argument("--condor", action='store_true', help="Flag which indicates we should submit jobs to Condor")

  # prepare data
  parser.add_argument("--data", help="CSV", required=True)
  parser.add_argument("--nodelist", help="Contains node identifiers for data and manifolds", required=True)

  # randomization
  # (run kegg.R then run randomize_network.R then prepare two separate files manifolds_file_true 
  # and manifolds_file_random which contain filepaths to some of the kegg and some of randomized kegg, respectively)
  # TODO trouble is the selection of a subset of networks into a file and that randomize_network.R runs on directories
  # could be solved by specifying selection with KEGG pathway identifiers and a directory rather than by a filepath

  # factorization
  parser.add_argument("--manifolds-file-true", required=True)
  parser.add_argument("--manifolds-file-random", required=True)
  parser.add_argument("--gamma", default="1.0")
  parser.add_argument("--k-latent", default="6")
  # other arguments used in factorization: outdir, nodelist, data

  # evaluation
  # (no additional arguments)
  args = parser.parse_args()

  job_graph = nx.DiGraph()
  job_id = 0

  # factorization
  # two branches: nmf_pathway on true pathways and nmf_pathway on randomized pathways
  # 1) nmf_pathway on true - {{
  nmf_pathway_true_outdir = os.path.join(args.outdir, "nmf_pathway_true")
  os.mkdir(nmf_pathway_true_outdir)
  attrs = {
    'exe': "nmf_pathway.py",
    'args': ["--data", args.data, '--k-latent', args.k_latent, "--manifolds-file", args.manifolds_file_true, "--nodelist", args.nodelist, "--node-attribute", "name", "--gamma", args.gamma, "--outdir", nmf_pathway_true_outdir],
    'out': os.path.join(nmf_pathway_true_outdir, "nmf_pathway.out"),
    'err': os.path.join(nmf_pathway_true_outdir, "nmf_pathway.err")
  }
  if(args.rng_seed is not None):
    attrs['args'] += ['--seed', args.rng_seed]
  prmf_true_gene_by_latent_fp = os.path.join(nmf_pathway_true_outdir, "V.csv")
  prmf_true_job_id = job_id
  job_graph.add_node(prmf_true_job_id, attrs)
  job_graph.add_edge(diffusion_job_id, job_id)
  job_id += 1

  # evaluation
  attrs = {
    'exe': "evaluate_screen_sim.py",
    'args': ["--gene-by-latent", prmf_true_gene_by_latent_fp, "--nodelist", nodelist_fp, "--true-seeds", chosen_seeds_fp],
    'out': os.path.join(nmf_pathway_true_outdir, "evaluate.out"),
    'err': os.path.join(nmf_pathway_true_outdir, "evaluate.err")
  }
  job_graph.add_node(job_id, attrs)
  job_graph.add_edge(job_id-1, job_id)
  job_id += 1
  # }} - nmf_pathway

  # 2) nmf_pathway on random - {{
  nmf_pathway_random_outdir = os.path.join(args.outdir, "nmf_pathway_random")
  os.mkdir(nmf_pathway_random_outdir)
  attrs = {
    'exe': "nmf_pathway.py",
    'args': ["--data", diffused_fp, '--k-latent', args.k_latent, "--manifolds-file", args.manifolds_file_random, "--nodelist", nodelist_fp, "--node-attribute", "name", "--gamma", args.gamma, "--outdir", nmf_pathway_random_outdir],
    'out': os.path.join(nmf_pathway_random_outdir, "nmf_pathway.out"),
    'err': os.path.join(nmf_pathway_random_outdir, "nmf_pathway.err")
  }
  if(args.rng_seed is not None):
    attrs['args'] += ['--seed', args.rng_seed]
  prmf_random_gene_by_latent_fp = os.path.join(nmf_pathway_random_outdir, "V.csv")
  prmf_random_job_id = job_id
  job_graph.add_node(prmf_random_job_id, attrs)
  job_graph.add_edge(diffusion_job_id, job_id)
  job_id += 1

  # evaluation
  attrs = {
    'exe': "evaluate_screen_sim.py",
    'args': ["--gene-by-latent", prmf_random_gene_by_latent_fp, "--nodelist", nodelist_fp, "--true-seeds", chosen_seeds_fp],
    'out': os.path.join(nmf_pathway_random_outdir, "evaluate.out"),
    'err': os.path.join(nmf_pathway_random_outdir, "evaluate.err")
  }
  job_graph.add_node(job_id, attrs)
  job_graph.add_edge(job_id-1, job_id)
  job_id += 1
  # }} - 

  # plot
  # TODO rework name of arguments
  plot_outdir = os.path.join(args.outdir, 'pr_curves')
  os.mkdir(plot_outdir)
  attrs = {
    'exe': 'plot_pr_curve.py',
    'args': [
      '--gene-by-latent-csvs', prmf_random_gene_by_latent_fp, prmf_true_gene_by_latent_fp, 
      '--labels', 'PRMF_random', 'PRMF_real',
      '--nodelist', nodelist_fp, '--true-seeds', chosen_seeds_fp, '--outdir', plot_outdir
    ],
    'out': os.path.join(plot_outdir, 'plot.out'),
    'err': os.path.join(plot_outdir, 'plot.err')
  }
  job_graph.add_node(job_id, attrs)
  # run after nmf and nmf_pathway
  job_graph.add_edge(prmf_true_job_id, job_id)
  job_graph.add_edge(prmf_random_job_id, job_id)
  job_id += 1

  condor = False
  if args.condor:
    condor = True
  job_ids = script_utils.run_digraph(args.outdir, job_graph, condor=condor)

if __name__ == "__main__":
  main()
