#!/usr/bin/env python
import argparse, sys
import os, os.path
import networkx as nx
from factorlib import script_utils

def main():
  parser = argparse.ArgumentParser(description=
"""
Evalute nmf_pathway.py by simulating gene lists and compare against nmf_init.py
""")
  parser.add_argument("--rng-seed", help="Seed for random number generators", default=None)

  # simulation
  parser.add_argument("--n-gene-lists", help="Number of gene lists to simulate", type=int, default=6)
  parser.add_argument("--nodelist", help="Universe of node identifiers and an ordering on those identifiers", required=True)
  parser.add_argument("--seed-lists", required=True, nargs='+') # TODO future versions manifolds should be distinct from node lists derived from pathways
  parser.add_argument("--outdir", required=True)
  parser.add_argument("--simulator")
  parser.add_argument("--noise-pr", default="0.05")

  # diffusion
  parser.add_argument("--network", required=True)
  # other arguments used in diffusion: nodelist

  # factorization
  parser.add_argument("--manifolds-file", required=True)
  parser.add_argument("--gamma", default="1.0")
  parser.add_argument("--k-latent", default="6")
  # other arguments used in factorization: outdir, nodelist, data

  # evaluation
  # (no additional arguments)
  args = parser.parse_args()

  job_graph = nx.DiGraph()
  job_id = 0

  # simulation
  attrs = {
    'exe': "simulate_screens.py",
    'args': ["--seed-lists"] + args.seed_lists + ["--n-gene-lists", str(args.n_gene_lists), "--nodelist", args.nodelist, "--outdir", args.outdir, '--simulator', args.simulator, '--noise-pr', args.noise_pr],
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
  diffusion_job_id = job_id
  job_graph.add_node(job_id, attrs)
  job_graph.add_edge(job_id-1, job_id)
  job_id += 1

  # factorization
  # two branches: nmf and nmf_pathway
  # 1) nmf - {{
  nmf_outdir = os.path.join(args.outdir, "nmf")
  os.mkdir(nmf_outdir)
  attrs = {
    'exe': "nmf.py",
    'args': ['--data', diffused_fp, '--k-latent', args.k_latent, '--outdir', nmf_outdir],
    'out': os.path.join(nmf_outdir, 'nmf.out'),
    'err': os.path.join(nmf_outdir, 'nmf.err')
  }
  if(args.rng_seed is not None):
    attrs['args'] += ['--seed', args.rng_seed]
  nmf_gene_by_latent_fp = os.path.join(nmf_outdir, "V.csv")
  nmf_job_id = job_id
  job_graph.add_node(nmf_job_id, attrs)
  job_graph.add_edge(diffusion_job_id, nmf_job_id)
  job_id += 1

  attrs = {
    'exe': "evaluate_simulation.py",
    'args': ["--gene-by-latent", nmf_gene_by_latent_fp, "--nodelist", args.nodelist, "--true-seeds", chosen_seeds_fp],
    'out': os.path.join(nmf_outdir, "evaluate.out"),
    'err': os.path.join(nmf_outdir, "evaluate.err")
  }
  job_graph.add_node(job_id, attrs)
  job_graph.add_edge(job_id-1, job_id)
  job_id += 1
  # }} - nmf

  # 2) nmf_pathway - {{
  nmf_pathway_outdir = os.path.join(args.outdir, "nmf_pathway")
  os.mkdir(nmf_pathway_outdir)
  attrs = {
    'exe': "nmf_pathway.py",
    'args': ["--data", diffused_fp, '--k-latent', args.k_latent, "--manifolds-file", args.manifolds_file, "--nodelist", args.nodelist, "--gamma", args.gamma, "--outdir", nmf_pathway_outdir],
    'out': os.path.join(nmf_pathway_outdir, "nmf_pathway.out"),
    'err': os.path.join(nmf_pathway_outdir, "nmf_pathway.err")
  }
  if(args.rng_seed is not None):
    attrs['args'] += ['--seed', args.rng_seed]
  prmf_gene_by_latent_fp = os.path.join(nmf_pathway_outdir, "V.csv")
  prmf_job_id = job_id
  job_graph.add_node(prmf_job_id, attrs)
  job_graph.add_edge(diffusion_job_id, job_id)
  job_id += 1

  # evaluation
  attrs = {
    'exe': "evaluate_simulation.py",
    'args': ["--gene-by-latent", prmf_gene_by_latent_fp, "--nodelist", args.nodelist, "--true-seeds", chosen_seeds_fp],
    'out': os.path.join(nmf_pathway_outdir, "evaluate.out"),
    'err': os.path.join(nmf_pathway_outdir, "evaluate.err")
  }
  job_graph.add_node(job_id, attrs)
  job_graph.add_edge(job_id-1, job_id)
  job_id += 1
  # }} - nmf_pathway


  # 3) PLIER - {{
  PLIER_outdir = os.path.join(args.outdir, "PLIER")
  os.mkdir(PLIER_outdir)
  attrs = {
    'exe': 'PLIER_wrapper.R',
    'args': ['--data', diffused_fp, '--nodelist', args.nodelist, '--k-latent', args.k_latent, '--pathways-file', args.manifolds_file, '--outdir', PLIER_outdir],
    'out': os.path.join(PLIER_outdir, "PLIER_wrapper.out"),
    'err': os.path.join(PLIER_outdir, "PLIER_wrapper.err")
  }
  if(args.rng_seed is not None):
    attrs['args'] += ['--seed', args.rng_seed]
  PLIER_gene_by_latent_fp = os.path.join(PLIER_outdir, "Z.csv")
  PLIER_job_id = job_id
  job_graph.add_node(PLIER_job_id, attrs)
  job_graph.add_edge(diffusion_job_id, PLIER_job_id)
  job_id += 1

  # evaluation
  attrs = {
    'exe': "evaluate_simulation.py",
    'args': ["--gene-by-latent", PLIER_gene_by_latent_fp, "--nodelist", args.nodelist, "--true-seeds", chosen_seeds_fp],
    'out': os.path.join(PLIER_outdir, "evaluate.out"),
    'err': os.path.join(PLIER_outdir, "evaluate.err")
  }
  job_graph.add_node(job_id, attrs)
  job_graph.add_edge(job_id-1, job_id)
  job_id += 1
  # }} - PLIER

  # plot
  plot_outdir = os.path.join(args.outdir, 'pr_curves')
  os.mkdir(plot_outdir)
  attrs = {
    'exe': 'plot_pr_curve.py',
    'args': ['--nmf-gene-by-latent', nmf_gene_by_latent_fp, '--prmf-gene-by-latent', prmf_gene_by_latent_fp, '--plier-gene-by-latent', PLIER_gene_by_latent_fp, '--nodelist', args.nodelist, '--true-seeds', chosen_seeds_fp, '--outdir', plot_outdir],
    'out': os.path.join(plot_outdir, 'plot.out'),
    'err': os.path.join(plot_outdir, 'plot.err')
  }
  job_graph.add_node(job_id, attrs)
  # run after nmf and nmf_pathway
  job_graph.add_edge(nmf_job_id, job_id)
  job_graph.add_edge(prmf_job_id, job_id)
  job_graph.add_edge(PLIER_job_id, job_id)
  job_id += 1

  condor = False
  job_ids = script_utils.run_digraph(args.outdir, job_graph, condor=condor)

if __name__ == "__main__":
  main()
