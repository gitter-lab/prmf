#!/usr/bin/env python
import argparse, sys
import os, os.path
import networkx as nx
from prmf import script_utils

def main():
  parser = argparse.ArgumentParser(description=
"""
Evalute nmf_pathway.py by simulating gene lists and compare against nmf_init.py
""")
  # run environment
  parser.add_argument("--rng-seed", help="Seed for random number generators", default=None)
  parser.add_argument("--condor", action='store_true')
  parser.add_argument("--dry-run", action='store_true')
  parser.add_argument("--do-nmf", default=True, type=script_utils.str2bool)
  parser.add_argument("--do-plier", default=True, type=script_utils.str2bool)
  parser.add_argument("--do-nbs", default=True, type=script_utils.str2bool)
  parser.add_argument("--do-cogaps", default=True, type=script_utils.str2bool)

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
  parser.add_argument("--cross-validation", help="Fraction of data in [0,1] to hold out and test reconstruction performance on")
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
    'err': os.path.join(args.outdir, "simulate_screens.err"),
    'env': 'prmf'
  }
  if(args.rng_seed is not None):
    attrs['args'] += ['--rng-seed', args.rng_seed]
  job_graph.add_node(job_id, attrs)
  simulation_job_id = job_id
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
    'err': os.path.join(args.outdir, "diffusion.err"),
    'env': 'prmf'
  }
  diffusion_job_id = job_id
  job_graph.add_node(job_id, attrs)
  job_graph.add_edge(job_id-1, job_id)
  job_id += 1

  # factorization
  # two branches: nmf and nmf_pathway
  # 1) nmf - {{
  if args.do_nmf:
    nmf_outdir = os.path.join(args.outdir, "nmf")
    script_utils.mkdir_p(nmf_outdir)
    attrs = {
      'exe': "nmf.py",
      'args': ['--data', diffused_fp, '--k-latent', args.k_latent, '--outdir', nmf_outdir],
      'out': os.path.join(nmf_outdir, 'nmf.out'),
      'err': os.path.join(nmf_outdir, 'nmf.err'),
      'env': 'prmf'
    }
    if(args.rng_seed is not None):
      attrs['args'] += ['--seed', args.rng_seed]
    if(args.cross_validation is not None):
      attrs['args'] += ['--cross-validation', args.cross_validation]
    nmf_gene_by_latent_fp = os.path.join(nmf_outdir, "V.csv")
    nmf_job_id = job_id
    job_graph.add_node(nmf_job_id, attrs)
    job_graph.add_edge(diffusion_job_id, nmf_job_id)
    job_id += 1

    attrs = {
      'exe': "evaluate_screen_sim.py",
      'args': ["--gene-by-latent", nmf_gene_by_latent_fp, "--nodelist", args.nodelist, "--true-seeds", chosen_seeds_fp],
      'out': os.path.join(nmf_outdir, "evaluate.out"),
      'err': os.path.join(nmf_outdir, "evaluate.err"),
      'env': 'prmf'
    }
    job_graph.add_node(job_id, attrs)
    job_graph.add_edge(job_id-1, job_id)
    job_id += 1
  # }} - nmf

  # 2) nmf_pathway - {{
  nmf_pathway_outdir = os.path.join(args.outdir, "nmf_pathway")
  script_utils.mkdir_p(nmf_pathway_outdir)
  attrs = {
    'exe': "nmf_pathway.py",
    'args': ["--data", diffused_fp, '--k-latent', args.k_latent, "--manifolds-file", args.manifolds_file, "--node-attribute", "name", "--nodelist", args.nodelist, "--gamma", args.gamma, "--outdir", nmf_pathway_outdir],
    'out': os.path.join(nmf_pathway_outdir, "nmf_pathway.out"),
    'err': os.path.join(nmf_pathway_outdir, "nmf_pathway.err"),
    'env': 'prmf'
  }
  if(args.rng_seed is not None):
    attrs['args'] += ['--seed', args.rng_seed]
  if(args.cross_validation is not None):
    attrs['args'] += ['--cross-validation', args.cross_validation]
  prmf_gene_by_latent_fp = os.path.join(nmf_pathway_outdir, "V.csv")
  prmf_job_id = job_id
  job_graph.add_node(prmf_job_id, attrs)
  job_graph.add_edge(diffusion_job_id, job_id)
  job_id += 1

  # evaluation
  attrs = {
    'exe': "evaluate_screen_sim.py",
    'args': ["--gene-by-latent", prmf_gene_by_latent_fp, "--nodelist", args.nodelist, "--true-seeds", chosen_seeds_fp],
    'out': os.path.join(nmf_pathway_outdir, "evaluate.out"),
    'err': os.path.join(nmf_pathway_outdir, "evaluate.err"),
    'env': 'prmf'
  }
  job_graph.add_node(job_id, attrs)
  job_graph.add_edge(job_id-1, job_id)
  job_id += 1
  # }} - nmf_pathway


  # 3) PLIER - {{
  if args.do_plier:
    PLIER_outdir = os.path.join(args.outdir, "PLIER")
    script_utils.mkdir_p(PLIER_outdir)
    attrs = {
      'exe': 'PLIER_wrapper.R',
      'args': ['--data', diffused_fp, '--nodelist', args.nodelist, '--k-latent', args.k_latent, '--pathways-file', args.manifolds_file, '--node-attribute', 'name', '--outdir', PLIER_outdir],
      'out': os.path.join(PLIER_outdir, "PLIER_wrapper.out"),
      'err': os.path.join(PLIER_outdir, "PLIER_wrapper.err"),
      'env': 'PLIER'
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
      'exe': "evaluate_screen_sim.py",
      'args': ["--gene-by-latent", PLIER_gene_by_latent_fp, "--nodelist", args.nodelist, "--true-seeds", chosen_seeds_fp],
      'out': os.path.join(PLIER_outdir, "evaluate.out"),
      'err': os.path.join(PLIER_outdir, "evaluate.err"),
      'env': 'prmf'
    }
    job_graph.add_node(job_id, attrs)
    job_graph.add_edge(job_id-1, job_id)
    job_id += 1
  # }} - PLIER

  # 4) NBS - {{
  # NOTE NBS does its own diffusion based on binary somatic mutation profiles so we pass the simulated hits rather than our diffused data
  if args.do_nbs:
    NBS_outdir = os.path.join(args.outdir, 'NBS')
    script_utils.mkdir_p(NBS_outdir)
    attrs = {
      'exe': 'pyNBS_wrapper.py',
      'args': ['--nodelist', args.nodelist, '--gene-lists'] + sim_list_fps + ['--network', args.network, '--k-latent', args.k_latent, '--outdir', NBS_outdir],
      'out': os.path.join(NBS_outdir, 'pyNBS_wrapper.out'),
      'err': os.path.join(NBS_outdir, 'pyNBS_wrapper.err'),
      'env': 'pyNBS'
    }
    NBS_job_id = job_id
    NBS_gene_by_latent_fp = os.path.join(NBS_outdir, "W.csv")
    job_graph.add_node(NBS_job_id, attrs)
    job_graph.add_edge(simulation_job_id, NBS_job_id)
    job_id += 1

    # evaluation
    attrs = {
      'exe': "evaluate_screen_sim.py",
      'args': ["--gene-by-latent", NBS_gene_by_latent_fp, "--nodelist", args.nodelist, "--true-seeds", chosen_seeds_fp],
      'out': os.path.join(NBS_outdir, "evaluate.out"),
      'err': os.path.join(NBS_outdir, "evaluate.err"),
      'env': 'prmf'
    }
    job_graph.add_node(job_id, attrs)
    job_graph.add_edge(NBS_job_id, job_id)
    job_id += 1
  # }} - NBS

  # 5) CoGAPS - {{
  if args.do_cogaps:
    CoGAPS_outdir = os.path.join(args.outdir, 'CoGAPS')
    script_utils.mkdir_p(CoGAPS_outdir)
    attrs = {
      'exe': 'CoGAPS_wrapper.R',
      'args': ['--data', diffused_fp, '--k-latent', args.k_latent, '--outdir', args.outdir],
      'out': os.path.join(CoGAPS_outdir, 'CoGAPS_wrapper.out'),
      'err': os.path.join(CoGAPS_outdir, 'CoGAPS_wrapper.err'),
      'env': 'CoGAPS'
    }
    CoGAPS_job_id = job_id
    CoGAPS_gene_by_latent_fp = os.path.join(CoGAPS_outdir, "P.csv")
    job_graph.add_node(CoGAPS_job_id, attrs)
    job_graph.add_edge(diffusion_job_id, CoGAPS_job_id)
    job_id += 1

    # evaluation
    attrs = {
      'exe': "evaluate_screen_sim.py",
      'args': ["--gene-by-latent", CoGAPS_gene_by_latent_fp, "--nodelist", args.nodelist, "--true-seeds", chosen_seeds_fp],
      'out': os.path.join(CoGAPS_outdir, "evaluate.out"),
      'err': os.path.join(CoGAPS_outdir, "evaluate.err"),
      'env': 'prmf'
    }
    job_graph.add_node(job_id, attrs)
    job_graph.add_edge(NBS_job_id, job_id)
    job_id += 1
  # }} - CoGAPS

  # plot
  plot_outdir = os.path.join(args.outdir, 'pr_curves')
  script_utils.mkdir_p(plot_outdir)

  gene_by_latent_csvs = []
  labels = []
  if args.do_nmf:
    labels.append('NMF')
    nmf_gene_by_latent_fp = os.path.join(nmf_outdir, "V.csv")
    gene_by_latent_csvs.append(nmf_gene_by_latent_fp)

  labels.append('PRMF')
  gene_by_latent_csvs.append(prmf_gene_by_latent_fp)

  if args.do_plier:
    labels.append("PLIER")
    PLIER_gene_by_latent_fp = os.path.join(PLIER_outdir, "Z.csv")
    gene_by_latent_csvs.append(PLIER_gene_by_latent_fp)

  if args.do_nbs:
    labels.append("NBS")
    NBS_gene_by_latent_fp = os.path.join(NBS_outdir, "W.csv")
    gene_by_latent_csvs.append(NBS_gene_by_latent_fp)

  if args.do_cogaps:
    labels.append("CoGAPS")
    CoGAPS_gene_by_latent_fp = os.path.join(CoGAPS_outdir, "P.csv")
    gene_by_latent_csvs.append(CoGAPS_gene_by_latent_fp)

  attrs = {
    'exe': 'plot_pr_curve.py',
    'args': [
      '--gene-by-latent-csvs'] + gene_by_latent_csvs +
      ['--labels'] + labels +
      ['--nodelist', args.nodelist, '--true-seeds', chosen_seeds_fp, '--outdir', plot_outdir],
    'out': os.path.join(plot_outdir, 'plot.out'),
    'err': os.path.join(plot_outdir, 'plot.err'),
    'env': 'prmf'
  }
  job_graph.add_node(job_id, attrs)
  # run after all methods
  if args.do_nmf:
    job_graph.add_edge(nmf_job_id, job_id)
  job_graph.add_edge(prmf_job_id, job_id)
  if args.do_plier:
    job_graph.add_edge(PLIER_job_id, job_id)
  if args.do_nbs:
    job_graph.add_edge(NBS_job_id, job_id)
  if args.do_cogaps:
    job_graph.add_edge(CoGAPS_job_id, job_id)
  job_id += 1

  condor = False
  if args.condor:
    condor = True
  job_ids = script_utils.run_digraph(args.outdir, job_graph, condor=condor, dry_run=args.dry_run)

if __name__ == "__main__":
  main()
