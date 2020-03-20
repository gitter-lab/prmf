#!/usr/bin/env python
import argparse, sys

def main():
  parser = argparse.ArgumentParser(description="""
Analyze PLIER output U matrix to determine which pathways are associated with the latent factors
""")
  parser.add_argument("--U", required=True)
  parser.add_argument("--outdir", required=True)
  args = parser.parse_args()

  # skip first row and col
  U = np.genfromtxt(fp, skip_header=1, delimiter=",")[:,1:]
  n_pathway, k_latent = U.shape
  k_to_pathway_and_scores = {}
  for k in range(k_latent):
    # there is an L1 penalty on U so many values should be 0, filter them out
    # and return list in order of top scoring pathways
    k_to_pathway_and_scores[k] = sorted(filter(lambda p: p[1] != 0, enumerate(U[:,k])), key=lambda p: p[1], reverse=True)

  for k in range(k_latent):
    pathway_and_scores = k_to_pathway_and_scores[k]

    print("{} -> {}".format(k, ','.join(map(

if __name__ == "__main__":
  main()
