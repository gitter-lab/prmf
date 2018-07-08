#!/usr/bin/env python
import argparse, sys
import os, os.path
import networkx as nx
from factorlib import script_utils

def main():
  parser = argparse.ArgumentParser(description=
"""
Evalute nmf_pathway.py by simulating gene lists
""")
  parser.add_argument("--rng-seed", help="Seed for random number generators", default=None)

  # simulation
  parser.add_argument("--n-gene-lists", help="Number of gene lists to simulate", type=int, default=6)
  parser.add_argument("--nodelist", help="Universe of node identifiers and an ordering on those identifiers", required=True)
  parser.add_argument("--seed-lists", required=True, nargs='+') # TODO future versions manifolds should be distinct from node lists derived from pathways
  parser.add_argument("--outdir", required=True)
  parser.add_argument("--simulator")

  # diffusion
  parser.add_argument("--network", required=True)
  #parser.add_argument("--nodelist")

  # factorization
  #parser.add_argument("--data", required=True)
  parser.add_argument("--manifolds", required=True, nargs='+')
  #parser.add_argument("--outdir")
  #parser.add_argument("--nodelist")

  # evaluation
  #
  args = parser.parse_args()

  job_graph = nx.DiGraph()
  job_id = 0

  # simulation
  attrs = {
    'exe': "simulate_screens.py",
    'args': ["--seed-lists"] + args.seed_lists + ["--n-gene-lists", str(args.n_gene_lists), "--nodelist", args.nodelist, "--outdir", args.outdir, '--simulator', args.simulator],
    'out': os.path.join(args.outdir, "simulate_screens.out"),
    'err': os.path.join(args.outdir, "simulate_screens.err")
  }
  if(args.rng_seed is not None):
    attrs['args'] += ['--rng-seed', args.rng_seed]
  job_graph.add_node(job_id, attrs)
  job_id += 1
  sim_list_fps = []
  for i in range(args.n_gene_lists):
    sim_list_fps.append(os.path.join(args.outdir, "sim_list_{}.txt".format(i+1)))
  chosen_seeds_fp = os.path.join(args.outdir, "chosen_seeds.txt")

  # diffusion
  diffused_fp = os.path.join(args.outdir, "diffused.csv")
  attrs = {
    'exe': "diffusion.py",
    'args': ["--network", args.network, "--nodelist", args.nodelist, "--gene-lists"] + sim_list_fps +  ["--diffused", diffused_fp],
    'out': os.path.join(args.outdir, "diffusion.out"),
    'err': os.path.join(args.outdir, "diffusion.err")
  }
  job_graph.add_node(job_id, attrs)
  job_graph.add_edge(job_id-1, job_id)
  job_id += 1

  # factorization
  attrs = {
    'exe': "nmf_pathway.py",
    'args': ["--data", diffused_fp, "--manifolds"] + args.manifolds + ["--nodelist", args.nodelist, "--outdir", args.outdir],
    'out': os.path.join(args.outdir, "nmf_pathway.out"),
    'err': os.path.join(args.outdir, "nmf_pathway.err")
  }
  if(args.rng_seed is not None):
    attrs['args'] += ['--seed', args.rng_seed]
  gene_by_latent_fp = os.path.join(args.outdir, "V.csv")
  job_graph.add_node(job_id, attrs)
  job_graph.add_edge(job_id-1, job_id)
  job_id += 1

  # evaluation
  # TODO new version
  attrs = {
    'exe': "evaluate_simulation.py",
    'args': ["--gene-by-latent", gene_by_latent_fp, "--nodelist", args.nodelist, "--true-seeds", chosen_seeds_fp],
    'out': os.path.join(args.outdir, "evaluate.out"),
    'err': os.path.join(args.outdir, "evaluate.err")
  }
  job_graph.add_node(job_id, attrs)
  job_graph.add_edge(job_id-1, job_id)
  job_id += 1

  condor = False
  job_ids = script_utils.run_digraph(args.outdir, job_graph, condor=condor)

if __name__ == "__main__":
  main()
