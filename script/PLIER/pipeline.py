#!/usr/bin/env python
import argparse, sys
import os, os.path
import networkx as nx
from factorlib import script_utils

if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument("--n-runs", default=1, type=int)
  parser.add_argument("--outdir")
  args = parser.parse_args()

  job_graph = nx.DiGraph()
  job_id = 0

  eval_job_ids = []
  for i in range(args.n_runs):
    # simluation
    run_outdir = os.path.join(args.outdir, "run{}".format(i))
    script_utils.mkdir_p(run_outdir)
    attrs = {
      'exe': 'simulate_PLIER.py',
      'args': ['-o', run_outdir, '--seed', str(i)],
      'out': os.path.join(run_outdir, 'generate.out'),
      'err': os.path.join(run_outdir, 'generate.err')
    }
    simulate_job_id = job_id
    job_graph.add_node(simulate_job_id, attrs)
    job_id += 1

    # PLIER
    PLIER_outdir = os.path.join(run_outdir, 'PLIER')
    script_utils.mkdir_p(PLIER_outdir)
    attrs = {
      'exe': 'PLIER_wrapper.R',
      'args': ['--data', os.path.join(run_outdir, 'Y.csv'), '--nodelist', os.path.join(run_outdir, 'nodelist.txt'), '--pathways-file', os.path.join(run_outdir, 'pathways_file.txt'), '--outdir', PLIER_outdir, '--k-latent', str(30), '--seed', str(i)],
      'out': os.path.join(PLIER_outdir, 'PLIER.out'),
      'err': os.path.join(PLIER_outdir, 'PLIER.err')
    }
    PLIER_job_id = job_id
    job_graph.add_node(PLIER_job_id, attrs)
    job_graph.add_edge(simulate_job_id, job_id)
    job_id += 1

    attrs = {
      'exe': 'evaluate_mse.py',
      'args': ['--true-z', os.path.join(run_outdir, "Z.csv"), '--pred-z', os.path.join(PLIER_outdir, 'Z.csv')],
      'out': os.path.join(PLIER_outdir, 'evaluate_mse.out'),
      'err': os.path.join(PLIER_outdir, 'evaluate_mse.err')
    }
    PLIER_eval_job_id = job_id
    eval_job_ids.append(PLIER_eval_job_id)
    job_graph.add_node(PLIER_eval_job_id, attrs)
    job_graph.add_edge(PLIER_job_id, PLIER_eval_job_id)
    job_id += 1

    # NMF
    NMF_outdir = os.path.join(run_outdir, 'NMF')
    script_utils.mkdir_p(NMF_outdir)
    attrs = {
      'exe': 'nmf.py',
      'args': ['--data', os.path.join(run_outdir, "Y.csv"), '--k-latent', str(30), '--seed', str(i), '--outdir', NMF_outdir],
      'out': os.path.join(NMF_outdir, 'nmf.out'),
      'err': os.path.join(NMF_outdir, 'nmf.err')
    }
    NMF_job_id = job_id
    job_graph.add_node(NMF_job_id, attrs)
    job_graph.add_edge(simulate_job_id, NMF_job_id)
    job_id += 1

    attrs = {
      'exe': 'evaluate_mse.py',
      'args': ['--true-z', os.path.join(run_outdir, 'Z.csv'), '--pred-z', os.path.join(NMF_outdir, 'V.csv')],
      'out': os.path.join(NMF_outdir, 'evaluate_mse.out'),
      'err': os.path.join(NMF_outdir, 'evaluate_mse.err')
    }
    NMF_eval_job_id = job_id
    eval_job_ids.append(NMF_eval_job_id)
    job_graph.add_node(NMF_eval_job_id, attrs)
    job_graph.add_edge(NMF_job_id, NMF_eval_job_id)
    job_id += 1

  attrs = {
    'exe': 'plot_mse.py',
    'args': ['-i', args.outdir, '-o', args.outdir],
    'out': os.path.join(args.outdir, 'plot_mse.out'),
    'err': os.path.join(args.outdir, 'plot_mse.err')
  }
  plot_job_id = job_id
  job_graph.add_node(plot_job_id, attrs)
  for eval_job_id in eval_job_ids:
    job_graph.add_edge(eval_job_id, plot_job_id)
  job_id += 1

  condor = False
  job_ids = script_utils.run_digraph(args.outdir, job_graph, exit_on_err=False, condor=condor)
