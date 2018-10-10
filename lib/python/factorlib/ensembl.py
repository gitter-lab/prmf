"""
Ensembl.org recommends the use of biomart to map HGNC, ENSG, ENST, ENSP
Currently these are implemented in the R script dl_ensembl_map
TODO add biomart python API functions here
"""
import sys
import re

def parse_mapping_tsv(fh, key_index=0, value_index=1, delim='\t'):
  """
  Parse a TSV file like the file at 
    ftp://ftp.ensembl.org/pub/release-94/tsv/homo_sapiens/Homo_sapiens.GRCh38.94.refseq.tsv.gz
  into a mapping 
  """
  rv = {}
  if type(fh) == 'str':
    fh = open(fh, 'r')
  for line in fh:
    line = line.rstrip()
    words = line.split(delim)
    key = words[key_index]
    value = None
    if(len(words) > value_index):
      value = words[value_index]
    else:
      # TODO warn?
      pass
    if value is not None:
      if key not in rv:
        rv[key] = set()
      rv[key].add(value)
  return rv

def map_hgnc_to_ensps(fh):
  """
  Return map from HGNC symbol to many ENSPs

  Parameters
  ----------
  fh : file-like 
    database file from dl_ensembl_map: 4 column TSV of HGNC, ENSG, ENST, ENSP

  Returns
  -------
  rv : dict
    mapping of HGNC symbol to a set of ENSPs
  """
  return parse_mapping_tsv(fh, key_index=0, value_index=3)

def filter_one_to_many(hgnc_to_ensps_map, G):
  """
  Return a map from HGNC symbol to one ENSP. The one is chosen to be the first ENSP in G.
  ENSPs are examined in sorted order.

  Parameters
  ----------
  hgnc_to_ensps_map : dict
    return value of map_hgnc_to_ensps

  copy_map : dict
    map hgnc to the number of ENSP in the map that also appear in G

  G : nx.Graph
    protein-protein interaction network such as STRING
  """
  rv = {}
  copy_map = {}
  for hgnc, ensps in hgnc_to_ensps_map.items():
    ensps = sorted(ensps)
    any_ensp_in_G = False
    for ensp in ensps:
      if ensp in G:
        any_ensp_in_G = True
        if hgnc in rv:
          # then this is the second ENSP in G
          if hgnc in copy_map:
            copy_map[hgnc] += 1
          else:
            copy_map[hgnc] = 2
        else:
          rv[hgnc] = ensp
    if not any_ensp_in_G:
      # then use the first ensp in sorted order
      if len(ensps) > 0:
        rv[hgnc] = ensps[0]

  return rv, copy_map

def map_ensp_to_hgnc(fh):
  rv = {}
  for line in fh:
    line = line.rstrip()
    words = line.split('\t')
    hgnc_symbol = words[0]
    ensp_id = None
    if(len(words) >= 4):
      ensp_id = words[3]
    if ensp_id is not None:
      if ensp_id in rv:
        sys.stderr.write("Encountered ENSP id {} more than once\n".format(ensp_id))
      rv[ensp_id] = hgnc_symbol
  return rv

def apply_map(word_map, gene_list_fh):
  """
  Read words (one per line) from gene_list_fh and apply word_map.
  Return mapped words as a set.

  Parameters
  ----------
  word_map : dict
    mapping of HGNC symbol to set of ENSPs

  gene_list_fh : file-like
    newline delimited file of HGNC symbols

  Returns
  -------
  rv : list
    sorted list of mapped ENSP identifiers
  """
  rv = set()
  for line in gene_list_fh:
    hgnc_symbol = line.rstrip()
    ensp_ids = word_map.get(hgnc_symbol)
    if ensp_ids is None:
      sys.stderr.write("Unrecognized HGNC symbol: {}\n".format(hgnc_symbol))
    else:
      rv = rv.union(ensp_ids)
  rv = sorted(rv)
  return rv

def apply_mapping(word_map, io_pairs):
  """
  For infile, outfile pairs in <io_pairs>, read words separated by PCRE non-word-character 
  regexp '\W', apply <word_map>, and write mapped words to outfile. Note many words in infile
  do not need to appear in the <word_map> and will be left unmapped.
  """
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
              ensp = word_map.get(word)
              if ensp is None:
                ofh.write(word + match_data.group(0))
              else:
                ofh.write(ensp + match_data.group(0))
              word = ""
            else:
              word += char
