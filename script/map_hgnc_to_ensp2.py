#!/usr/bin/env python
import os, os.path
import argparse, sys
import re
from factorlib import ensembl
from factorlib.string_db import parse_string_fh
from factorlib.script_utils import add_file_and_dir_args
from factorlib.script_utils import check_file_and_dir_args

def add_file_and_dir_args(parser):
  parser.add_argument("--infiles", nargs='+', type=str)
  parser.add_argument("--outfiles", nargs='+', type=str)
  parser.add_argument("--indir", "-i", type=str)
  parser.add_argument("--outdir", "-o", type=str)

def check_file_and_dir_args(args):
  # require file mode or directory mode
  do_file = False
  do_dir = False
  if args.infiles is not None and args.outfiles is not None:
    do_file = True
    if len(args.infiles) != len(args.outfiles):
      sys.stderr.write("Must provide the same number of infiles and outfiles\n")
      sys.exit(23)
  if not do_file and (args.indir is not None and args.outdir is not None):
    do_dir = True
  if not do_file and not do_dir:
    sys.stderr.write("Must provide --infiles and --outfiles or provide --indir and --outdir\n")
    sys.exit(22)
  io_pairs = []
  if do_dir:
    for ifn in os.listdir(args.indir):
      ifp = os.path.join(args.indir, ifn)
      ofp = os.path.join(args.outdir, ifn)
      io_pairs.append((ifp, ofp))
  else:
    # then do_file
    io_pairs = list(zip(args.infiles, args.outfiles))
  return io_pairs

def apply_mapping(word_map, io_pairs):
  sep_re = re.compile('\W')
  word = ""
  for infile, outfile in io_pairs:
    with open(infile, "r") as ifh:
      with open(outfile, "w") as ofh:
        for line in ifh:
          for char in line:
            match_data = sep_re.match(char)
            if(match_data is not None):
              # then process word
              ensp = hgnc_to_ensp_map.get(word)
              if ensp is None:
                ofh.write(word + match_data.group(0))
              else:
                ofh.write(ensp + match_data.group(0))
              word = ""
            else:
              word += char

def main():
  parser = argparse.ArgumentParser(description="""
Translate HGNC words in <infile> to ENSP. If HGNC cannot be mapped, leave it.
""")
  parser.add_argument("--mapping-file", "-m", help="See dl_ensembl_map.R", required=True, type=argparse.FileType('r'))
  parser.add_argument("--network-file", "-n", help="ENSP protein interaction network e.g. STRING (many genes have evidence reported only for one of the proteins they encode)", required=True, type=argparse.FileType('r'))
  add_file_and_dir_args(parser)
  args = parser.parse_args()

  io_pairs = check_file_and_dir_args(args)

  hgnc_to_ensps_map = ensembl.map_hgnc_to_ensps(args.mapping_file)
  G = parse_string_fh(args.network_file)
  hgnc_to_ensp_map, copy_map = ensembl.filter_one_to_many(hgnc_to_ensps_map, G)
  for k, v in copy_map.iteritems():
    sys.stderr.write("[warning] {} ENSPs associated to HGNC:{} found in <network-file>\n".format(v, k))

  ensembl.apply_mapping(hgnc_to_ensp_map, io_pairs)

if __name__ == "__main__":
  main()
