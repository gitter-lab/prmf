#!/usr/bin/env python
import os, os.path
import argparse, sys
from prmf import ensembl
from prmf.string_db import parse_string_fh
from prmf.script_utils import add_file_and_dir_args, check_file_and_dir_args

def main():
  parser = argparse.ArgumentParser()
  parser.add_argument("--mapping-file", '-m', required=True, type=argparse.FileType('r'))
  parser.add_argument("--network-file", "-n", help="ENSP protein interaction network e.g. STRING (many genes have evidence reported only for one of the proteins they encode)", required=True, type=argparse.FileType('r'))
  add_file_and_dir_args(parser)
  args = parser.parse_args()

  io_pairs = check_file_and_dir_args(args)

  hgnc_to_ensps_map = ensembl.parse_mapping_tsv(args.mapping_file, key_index=0, value_index=2)
  G = parse_string_fh(args.network_file)
  ensg_to_ensp_map, copy_map = ensembl.filter_one_to_many(hgnc_to_ensps_map, G)
  for k, v in copy_map.items():
    sys.stderr.write("[warning] {} ENSPs associated to {} found in <network-file>\n".format(v, k))

  ensembl.apply_mapping(ensg_to_ensp_map, io_pairs)

if __name__ == "__main__":
  main()
