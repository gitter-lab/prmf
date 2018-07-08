#!/usr/bin/env python
import argparse
import os.path
import numpy as np
import factorlib as fl
from factorlib import ensembl
from factorlib.string_db import parse_string_fh

def main():
  parser = argparse.ArgumentParser(description="""
Prepare Achilles data from http://portals.broadinstitute.org/achilles/datasets/18/download
for input to nmf_pathway:
  (1) Map Achilles HUGO/HGNC gene symbols to ENSP for compatibility with global protein interaction network from STRINGdb
  (2) Construct a master nodelist by aggregating all node identifier present in Achilles, pathways, and STRINGdb
""")
  parser.add_argument("--gene-dependency", help="CSV file provided by Achilles")
  parser.add_argument("--mapping-file")
  parser.add_argument("--network-file")
  parser.add_argument("--pathways-dir")
  parser.add_argument("--outdir")
  args = parser.parse_args()

  hgnc_to_ensps_map = ensembl.map_hgnc_to_ensps(args.mapping_file)
  G = parse_string_fh(open(args.network_file))
  hgnc_to_ensp_map, copy_map = ensembl.filter_one_to_many(hgnc_to_ensps_map, G)

  csv_fp = os.path.join(args.outdir, "avana_data.csv")
  nodelist_fp = os.path.join(args.outdir, "avana_nodelist.csv")

  Gs = fl.parse_pathways(args.pathways_dir)
  data, row_names, col_names = fl.parse_achilles(args.gene_dependency)

  # map column names
  col_names_ensp_or_hgnc = []
  for hgnc in col_names:
    ensp = hgnc_to_ensp_map.get(hgnc)
    if ensp is None:
      col_names_ensp_or_hgnc.append(hgnc)
    else:
      col_names_ensp_or_hgnc.append(ensp)

  lists = list(map(lambda H: H.nodes(), Gs))
  lists.append(col_names_ensp_or_hgnc)
  lists.append(G.nodes())
  nodelist = fl.prepare_nodelist(lists)

  # embed data in len(nodelist)-dimensional space
  data_embeded = fl.embed_arr(nodelist, col_names, data)

  with open(nodelist_fp, 'w') as fh:
    fh.write('\n'.join(nodelist))
  np.savetxt(csv_fp, data_embeded, delimiter=",")

if __name__ == "__main__":
  main()
