#!/usr/bin/env python
import sys, argparse
import os, os.path
import networkx as nx
import factorlib as fl
from factorlib import script_utils

def main():
  parser = argparse.ArgumentParser(description="""
Generative simulation pipeline to benchmark PRMF against NMF, PLIER, CoGAPS, NBS
""")
  parser.add_argument("--n-simulations", "-n", default=2, type=int, help='Number of simulations to run: default 2')
  parser.add_argument("--outdir", help='Directory to write results to', required=True)
  parser.add_argument("--seed", help="Seed for random number generators", default=None)
  parser.add_argument("--condor", action='store_true', help="Run the pipeline on HTCondor")
  parser.add_argument("--dry-run", action='store_true', help="Report which commands will be run but don't actually run then")
  parser.add_argument("--do-nmf", default=True, help="If true, include NMF in the benchmarking; default true", type=script_utils.str2bool)
  parser.add_argument("--do-prmf", default=True, help="If true, include PRMF in the benchmarking; default true", type=script_utils.str2bool)
  parser.add_argument("--do-plier", default=True, help="If true, include PLIER in the benchmarking; default true", type=script_utils.str2bool)
  args = parser.parse_args()
  fl.script_utils.log_script(sys.argv)

  job_graph = nx.DiGraph()
  job_id = 0

  # TODO update use of seed
  nmf_job_ids = []
  prmf_job_ids = []
  plier_job_ids = []
  for i in range(args.n_simulations):
    outdir = os.path.join(args.outdir, 'sim{}'.format(i))
    script_utils.mkdir_p(outdir) 
    attrs = {
      'exe': 'nmf_pathway_sim_gen.py',
      'args': ['--outdir', outdir],
      'out': os.path.join(outdir, 'nmf_pathway_sim_gen.out'),
      'err': os.path.join(outdir, 'nmf_pathway_sim_gen.err'),
      'env': 'prmf'
    }
    if args.seed is not None:
      attrs['args'] += ['--seed', args.seed]
    data = os.path.join(outdir, 'X.csv')
    sim_sample_by_latent = os.path.join(outdir, 'U.csv')
    sim_gene_by_latent = os.path.join(outdir, 'V.csv')
    pathways = list(map(lambda x: os.path.join(outdir, 'pathway{}.graphml'.format(x)), range(1000)))
    pathways_file = os.path.join(outdir, 'pathways_file.txt')

    job_graph.add_node(job_id, attrs)
    sim_job_id = job_id
    job_id += 1

    nmf_job_id = None
    if args.do_nmf:
      nmf_outdir = os.path.join(outdir, 'nmf')
      script_utils.mkdir_p(nmf_outdir)

      # NMF
      attrs = {
        'exe': 'nmf.py',
        'args': ['--data', data, '--outdir', nmf_outdir, '--k-latent', '30'],
        'out': os.path.join(nmf_outdir, 'nmf.out'),
        'err': os.path.join(nmf_outdir, 'nmf.err'),
        'env': 'prmf'
      }
      if args.seed is not None:
        attrs['args'] += ['--seed', args.seed]
      nmf_gene_by_latent = os.path.join(nmf_outdir, 'V.csv')
      job_graph.add_node(job_id, attrs)
      job_graph.add_edge(sim_job_id, job_id)
      nmf_job_id = job_id
      nmf_job_ids.append(nmf_job_id)
      job_id += 1 

    prmf_job_id = None
    if args.do_prmf:
      prmf_outdir = os.path.join(outdir, 'prmf')
      script_utils.mkdir_p(prmf_outdir)

      # PRMF
      attrs = {
        'exe': 'prmf.py',
        'args': ['--data', data, '--manifolds'] + pathways + ['--node-attribute', 'name', '--k-latent', '30', '--outdir', prmf_outdir],
        'out': os.path.join(prmf_outdir, 'prmf.out'),
        'err': os.path.join(prmf_outdir, 'prmf.err'),
        'env': 'prmf'
      }
      job_graph.add_node(job_id, attrs)
      job_graph.add_edge(sim_job_id, job_id)
      prmf_job_id = job_id
      prmf_job_ids.append(prmf_job_id)
      job_id += 1

    plier_job_id = None
    if args.do_plier:
      plier_outdir = os.path.join(outdir, 'plier')
      script_utils.mkdir_p(plier_outdir)
      attrs = {
        'exe': 'PLIER_wrapper.R',
        'args': ['--data', data, '--pathways-file', pathways_file, '--k-latent', '30', '--node-attribute', 'name', '--L1', '50', '--L2', '50', '--outdir', plier_outdir],
        'out': os.path.join(plier_outdir, 'PLIER_wrapper.out'),
        'err': os.path.join(plier_outdir, 'PLIER_wrapper.err'),
        'env': 'prmf'
      }
      job_graph.add_node(job_id, attrs)
      job_graph.add_edge(sim_job_id, job_id)
      plier_job_id = job_id
      plier_job_ids.append(plier_job_id)
      job_id += 1

  # evaluation
  eval_outdir = os.path.join(args.outdir, 'eval')
  script_utils.mkdir_p(eval_outdir)
  attrs = {
    'exe': 'gen_sim_eval.py',
    'args': ['--indir', args.outdir, '--outdir', eval_outdir],
    'out': os.path.join(eval_outdir, 'gen_sim_eval.out'),
    'err': os.path.join(eval_outdir, 'gen_sim_eval.err'),
    'env': 'prmf'
  }
  job_graph.add_node(job_id, attrs)
  eval_job_id = job_id
  for nmf_job_id in nmf_job_ids:
    job_graph.add_edge(nmf_job_id, eval_job_id)
  for prmf_job_id in prmf_job_ids:
    job_graph.add_edge(prmf_job_id, eval_job_id)
  for plier_job_id in plier_job_ids:
    job_graph.add_edge(plier_job_id, eval_job_id)
  job_id += 1

  condor = False
  if args.condor:
    condor = True
  job_ids = script_utils.run_digraph(args.outdir, job_graph, condor=condor, dry_run=args.dry_run)

if __name__ == "__main__":
  main()
