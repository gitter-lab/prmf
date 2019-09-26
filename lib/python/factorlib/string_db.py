"""
Functions to parse STRING protein-protein interaction database:
head -n2 9606.protein.links.full.v10.5.txt
protein1 protein2 neighborhood neighborhood_transferred fusion cooccurence homology coexpression coexpression_transferred experiments experiments_transferred database database_transferred textmining textmining_transferred combined_score
9606.ENSP00000000233 9606.ENSP00000263431 0 0 0 0 0 0 53 0 176 0 0 0 128 260
"""
import networkx as nx

def parse_line_abc(line, score_index=-1):
  """
  Filter STRING data to just protein,protein,confidence and remove taxonomy code from protein identifiers
  """
  line = line.rstrip()
  words = line.split()
  p1 = words[0]
  p2 = words[1]
  # remove human taxonomy code prefix
  p1 = trim_taxonomy_code(p1)
  p2 = trim_taxonomy_code(p2)
  score = words[score_index]
  conf = int(score) / 1000.0
  return (p1, p2, conf)

def trim_taxonomy_code(protein_id):
  return protein_id[5:]

def parse_string_fh(fh, threshold=0.0, score='combined_score'):
  """
  Parameters
  ----------
  fh : file-like
    STRING db database file

  threshold : float
    edge weight confidence threshold expressed as a value in [0,1]
    include all edges with a confidence greater than or equal to the threshold

  score : str
    the name of a column in the string format to use as the edge confidence score
    default: combined_score
    Can be "neighborhood" "fusion" "cooccurence" "coexpression" "experimental" "database" 
    "textmining" or "combined_score"

  Returns
  -------
  G : nx.Graph
    protein-protein interaction network with Ensembl protein ids as node ids
  """
  # skip header
  line_no = 1
  header = fh.__next__()
  columns = header.strip().split()
  if score not in columns[2:]:
    # TODO exception class
    raise Exception("Invalid score column {}; must be one of {}".format(score, ", ".join(columns[2:])))
  score_index = columns.index(score)

  G = nx.Graph()
  for line in fh:
    line_no += 1
    p1, p2, conf = parse_line_abc(line, score_index=score_index)
    if(conf >= threshold):
      G.add_edge(p1, p2, {'weight': conf})
  return G

def string_to_abc(ifh, ofh):
  """
  Convert string database to "abc"-format

  Parameters
  ----------
  ifh : file-like
    String database

  ofh : file-like
    Output file handle to write abc graph to
  """
  # skip header
  ifh.next()
  for line in ifh:
    p1, p2, conf = parse_line_abc(line)
    ofh.write('\t'.join(map(str, [p1, p2, conf])) + '\n')
