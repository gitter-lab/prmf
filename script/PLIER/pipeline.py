#!/usr/bin/env python
import argparse, sys
import os, os.path
from factorlib import script_utils

if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument("--outdir")
  args = parser.parse_args()

  job_graph = nx.DiGraph()
  job_id = 0

  # simluation
  attrs = {
    'exe': 'simulate_plier.py',
    'args': ['-o', args.outdir],
    'out': os.path.join(args.outdir, 'generate.out'),
    'err': os.path.join(args.outdir, 'generate.err')
  }
  simulate_job_id = job_id
  job_graph.add_node(simulate_job_id, attrs)
  job_id += 1

  # PLIER
  PLIER_outdir = os.path.join(args.outdir, 'PLIER')
  os.mkdir(PLIER_outdir)
  attrs = {
    'exe': 'PLIER_wrapper.R',
    'args': ['--data', os.path.join(args.outdir, 'Y.csv'), '--nodelist', os.path.join(args.outdir, 'nodelist.txt'), '--pathways-file', os.path.join(args.outdir, 'pathways_file.txt'), '--outdir', PLIER_outdir, '--k-latent', str(4), '--seed', str(1)],
    'out': os.path.join(PLIER_outdir, 'PLIER.out'),
    'err': os.path.join(PLIER_outdir, 'PLIER.err')
  }
  PLIER_job_id = job_id
  job_graph.add_node(PLIER_job_id, attrs)
  job_graph.add_edge(simulate_job_id, job_id)
  job_id += 1

  attrs = {
    'exe': 'evaluate_mse.py',
    'args': ['--true-z', os.path.join(args.outdir, "Z.csv"), '--pred-z', os.path.join(PLIER_outdir, 'Z.csv')],
    'out': os.path.join(PLIER_outdir, 'evaluate_mse.out'),
    'err': os.path.join(PLIER_outdir, 'evaluate_mse.err')
  }
  PLIER_eval_job_id = job_id
  job_graph.add_node(PLIER_eval_job_id, attrs)
  job_graph.add_edge(PLIER_job_id, PLIER_eval_job_id)
  job_id += 1

  # NMF
  NMF_outdir = os.path.join(args.outdir, 'NMF')
  os.mkdir(NMF_outdir)
  attrs = {
    'exe': 'nmf.py',
    'args': ['--data', os.path.join(args.outdir, "Y.csv"), '--k-latent', str(4), '--seed', str(1), '--outdir', NMF_outdir],
    'out': os.path.join(NMF_outdir, 'nmf.out'),
    'err': os.path.join(NMF_outdir, 'nmf.err')
  }
  NMF_job_id = job_id
  job_graph.add_node(NMF_job_id, attrs)
  job_graph.add_edge(simulate_job_id, NMF_job_id)
  job_id += 1

  attrs = {
    'exe': 'evaluate_mse.py',
    'args': ['--true-z', os.path.join(args.outdir, 'Z.csv'), '--pred-z', os.path.join(NMF_outdir, 'V.csv')],
    'out': os.path.join(NMF_outdir, 'evaluate_mse.out'),
    'err': os.path.join(NMF_outdir, 'evaluate_mse.err')
  }
  NMF_eval_job_id = job_id
  job_graph.add_node(NMF_eval_job_id, attrs)
  job_graph.add_edge(NMF_job_id, NMF_eval_job_id)
  job_id += 1

  condor = False
  job_ids = script_utils.run_digraph(args.outdir, job_graph, condor=condor)
