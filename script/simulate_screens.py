#!/usr/bin/env python
import sys, argparse
import os, os.path
import factorlib as fl
import numpy.random as nprand

# TODO rename pathways
def sample_pathways(args, nodelist):
  """
  Returns
  -------
  seed_lists_sample : list of 
  """
  subset_inds = nprand.choice(len(args.seed_lists), size=args.subset_sample_size)
  seed_lists_sample = []
  chosen_seed_fps = []
  for subset_ind in subset_inds:
    seed_list_fp = args.seed_lists[subset_ind]
    seed_list = fl.parse_seedlist(seed_list_fp)
    seed_lists_sample.append(seed_list)
    chosen_seed_fps.append(seed_list_fp)
  seed_lists_sample.append(nodelist) # final seed list is background <nodelist>

  seed_list_sizes = []
  for seed_list in seed_lists_sample:
    seed_list_sizes.append(len(seed_list))

  return seed_lists_sample, seed_list_sizes, chosen_seed_fps

def write_lists(args, gene_lists):
  for i, gene_list in enumerate(gene_lists):
    gene_list_fp = os.path.join(args.outdir, "sim_list_{}.txt".format(i+1))
    with open(gene_list_fp, 'w') as ofh:
      ofh.write("\n".join(gene_list))

def write_seeds(args, chosen_seed_fps):
  chosen_out_fp = os.path.join(args.outdir, "chosen_seeds.txt")
  with open(chosen_out_fp, 'w') as ofh:
    ofh.write("\n".join(chosen_seed_fps))

def simulate_whole(args):
  """
  In contrast to simulate_mixture, do not combine nodes from different pathways into a single gene list
  """
  nodelist = fl.parse_nodelist(args.nodelist)

  # first, sample subset of seed-lists
  seed_lists_sample, seed_list_sizes, chosen_seed_fps = sample_pathways(args, nodelist)

  # then, from the sampled seed list and the background <nodelist>, sample gene lists
  gene_lists = []
  for i in range(args.n_gene_lists):
    # sample from multinomial to determine number of elements coming from each seed list
    sample_size = nprand.binomial(args.gene_list_size, 1 - args.noise_pr)

    # then sample from seed list uniformly at random
    gene_list = set()
    seed_list = seed_lists_sample[i]
    if sample_size < len(seed_list):
      seed_list_inds = nprand.choice(len(seed_list), size=sample_size, replace=False)
      for seed_list_ind in seed_list_inds:
        gene_list.add(seed_list[seed_list_ind])
    else:
      # TODO warn?
      for seed in seed_list:
        gene_list.add(seed)

    # sample remaining from background
    # TODO catch error if size > len(nodelist)
    nodelist_inds = nprand.choice(len(nodelist), size=(args.gene_list_size - sample_size), replace=False)
    for nodelist_ind in nodelist_inds:
      gene_list.add(nodelist[nodelist_ind])

    gene_lists.append(sorted(gene_list))

  # write gene lists to file
  write_lists(args, gene_lists)
  write_seeds(args, chosen_seed_fps)

def simulate_mixture(args):
  nodelist = fl.parse_nodelist(args.nodelist)

  # first, sample subset of seed-lists
  seed_lists_sample, seed_list_sizes, chosen_seed_fps = sample_pathways(args, nodelist)

  # define multinomial parameters
  prs = [(1 - args.noise_pr)/(len(seed_lists_sample) - 1)] * (len(seed_lists_sample) - 1)
  prs.append(args.noise_pr)

  # then, from the sampled seed list and the background <nodelist>, sample gene lists
  gene_lists = []
  for i in range(args.n_gene_lists):
    # sample from multinomial to determine number of elements coming from each seed list
    sample_sizes = nprand.multinomial(args.gene_list_size, prs)

    # if we try to overdraw from a seed list, sample the remaining from background
    for j in range(len(sample_sizes)-1):
      seed_list_size = seed_list_sizes[j]
      sample_size = sample_sizes[j]
      diff = sample_size - seed_list_size
      if diff > 0:
        sample_sizes[j] -= diff
        sample_sizes[-1] += diff

    # then sample from seed lists uniformly at random
    gene_list = set()
    for j in range(len(seed_lists_sample)):
      seed_list = seed_lists_sample[j]
      sample_size = sample_sizes[j]
      seed_list_inds = nprand.choice(len(seed_list), size=sample_size)
      for seed_list_ind in seed_list_inds:
        gene_list.add(seed_list[seed_list_ind])

    gene_lists.append(sorted(gene_list))

  # write gene lists to file
  write_lists(args, gene_lists)
  write_seeds(args, chosen_seed_fps)

def main():
  parser = argparse.ArgumentParser(description="""
Sample a subset of seed lists which contains node identifiers from <nodelist> and use the
resulting node sets to construct the simulated gene lists.
""")
  parser.add_argument("--simulator", type=str, default='mixture', help="One of \"mixture\", \"whole\"")
  parser.add_argument("--seed-lists", type=str, nargs='+', required=True)
  parser.add_argument("--n-gene-lists", default=6, type=int, help="Number of gene lists to simulate")
  parser.add_argument("--nodelist", type=argparse.FileType('r'), required=True)
  parser.add_argument("--outdir", required=True, type=str)
  parser.add_argument("--subset-sample-size", "-p", type=int, default=5, help="Number of seed lists to sample")
  parser.add_argument("--gene-list-size", type=int, default=80)
  parser.add_argument("--rng-seed", default=None)
  parser.add_argument("--noise-pr", type=float, default=0.05, help="Probability of selecting a gene from the background gene list rather than from one of the selected pathways. Default = 0.05")
  args = parser.parse_args()

  # seed random number generator, dont use fixed seed unless provided
  rng_seed = args.rng_seed
  if(rng_seed is not None):
    rng_seed = int(rng_seed)
  nprand.seed(rng_seed)

  if args.simulator == 'mixture':
    simulate_mixture(args)
  elif args.simulator == 'whole':
    # TODO word for "not mixture"?
    simulate_whole(args)
  else:
    sys.stderr.write("Unrecognized simulator: \"{}\"\n".format(args.simulator))
    sys.exit(1)
    
if __name__ == "__main__":
  main()
