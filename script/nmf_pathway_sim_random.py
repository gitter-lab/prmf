#!/usr/bin/env python
import argparse, sys
import os, os.path
import networkx as nx
from prmf import script_utils

def main():
  parser = argparse.ArgumentParser(description=
"""
Evalute nmf_pathway.py against true pathways and randomized pathways
""")
  parser.add_argument("--rng-seed", help="Seed for random number generators", default=None)
  parser.add_argument("--condor", action='store_true', help="Flag which indicates we should submit jobs to Condor")

  # prepare networks
  parser.add_argument('--stringdb', help="STRINGdb database file e.g. 9606.protein.links.detailed.v10.txt", required=True)

  # simulation
  parser.add_argument("--n-gene-lists", help="Number of gene lists to simulate", type=int, default=6)
  parser.add_argument("--seed-lists", required=True, nargs='+') # TODO future versions manifolds should be distinct from node lists derived from pathways
  parser.add_argument("--outdir", required=True)
  parser.add_argument("--simulator", default='mixture') # TODO check option value
  parser.add_argument("--noise-pr", default="0.05")

  # randomization
  # (run kegg.R then run randomize_network.R then prepare two separate files manifolds_file_true 
  # and manifolds_file_random which contain filepaths to some of the kegg and some of randomized kegg, respectively)
  # TODO trouble is the selection of a subset of networks into a file and that kegg.R runs on directories
  # could be solved by specifying selection with KEGG pathway identifiers and a directory rather than by a filepath

  # diffusion
  # (no additional arguments)
  # other arguments used in diffusion: nodelist

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

  # prepare networks
  nodelist_fp = os.path.join(args.outdir, 'nodelist.txt')
  string_kegg_union_fp = os.path.join(args.outdir, 'string_kegg_union.graphml')
  manifold_fps = []
  with open(args.manifolds_file_true) as fh:
    for line in fh:
      line = line.rstrip()
      manifold_fps.append(line)
  attrs = {
    'exe': 'prepare_nodelist.py',
    'args': ['--stringdb', args.stringdb, '--graphmls'] +  manifold_fps + ['--out-nodelist', nodelist_fp, '--out-graph', string_kegg_union_fp, '--node-attribute', 'name'],
    'out': os.path.join(args.outdir, 'prepare_nodelist.out'),
    'err': os.path.join(args.outdir, 'prepare_nodelist.err')
  }
  job_graph.add_node(job_id, attrs)
  job_id += 1

  # simulation
  attrs = {
    'exe': "simulate_screens.py",
    'args': ["--seed-lists"] + args.seed_lists + ["--n-gene-lists", str(args.n_gene_lists), "--nodelist", nodelist_fp, "--outdir", args.outdir, '--simulator', args.simulator, '--noise-pr', args.noise_pr],
    'out': os.path.join(args.outdir, "simulate_screens.out"),
    'err': os.path.join(args.outdir, "simulate_screens.err")
  }
  if(args.rng_seed is not None):
    attrs['args'] += ['--rng-seed', args.rng_seed]
  job_graph.add_node(job_id, attrs)
  job_graph.add_edge(job_id-1, job_id)
  job_id += 1
  sim_list_fps = []
  for i in range(args.n_gene_lists):
    sim_list_fps.append(os.path.join(args.outdir, "sim_list_{}.txt".format(i+1)))
  chosen_seeds_fp = os.path.join(args.outdir, "chosen_seeds.txt")

  # diffusion
  diffused_fp = os.path.join(args.outdir, "diffused.csv")
  attrs = {
    'exe': "diffusion.py",
    'args': ["--network", string_kegg_union_fp, "--nodelist", nodelist_fp, "--gene-lists"] + sim_list_fps +  ["--diffused", diffused_fp],
    'out': os.path.join(args.outdir, "diffusion.out"),
    'err': os.path.join(args.outdir, "diffusion.err")
  }
  diffusion_job_id = job_id
  job_graph.add_node(job_id, attrs)
  job_graph.add_edge(job_id-1, job_id)
  job_id += 1

  # factorization
  # two branches: nmf_pathway on true pathways and nmf_pathway on randomized pathways
  # 1) nmf_pathway on true - {{
  nmf_pathway_true_outdir = os.path.join(args.outdir, "nmf_pathway_true")
  os.mkdir(nmf_pathway_true_outdir)
  attrs = {
    'exe': "nmf_pathway.py",
    'args': ["--data", diffused_fp, '--k-latent', args.k_latent, "--manifolds-file", args.manifolds_file_true, "--nodelist", nodelist_fp, "--node-attribute", "name", "--gamma", args.gamma, "--outdir", nmf_pathway_true_outdir],
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
      '--labels', 'PRMF_random', 'PRMF_true',
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
