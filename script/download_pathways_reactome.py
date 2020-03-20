#!/usr/bin/env python
import time
import re
import sys, argparse
import os, os.path
from urllib.parse import urlencode
import requests
import prmf

KEGG_URI_BASE = "http://identifiers.org/kegg.pathway"
REACTOME_URI_BASE = "http://identifiers.org/reactome/"
PC_GET_URI = "http://www.pathwaycommons.org/pc2/get"
PC_KEGG_URI = "http://pathwaycommons.org/pc2/kegg"
PC_TOP_PATHWAYS_URI = "http://www.pathwaycommons.org/pc2/top_pathways"

def reactome_uri(reactome_id):
  return REACTOME_URI_BASE + reactome_id

def kegg_uri(kegg_id):
  return KEGG_URI_BASE + "/" + kegg_id

def get_kegg_pathway_ids():
  """
  Return dict mapping pathway identifier to pathway title
  """
  rv = {}
  resp = requests.get('http://rest.kegg.jp/list/pathway/hsa')
  lines = resp.text.split('\n')
  # path:hsa00010 <TAB> Glycolysis / Gluconeogenesis - Homo sapiens (human)
  line_regexp = re.compile(r'path:(\w+)\W+(\w+)')
  for line in lines:
    match_data = line_regexp.match(line)
    if match_data is not None:
      rv[match_data.group(1)] = match_data.group(2)
  return rv

def get_reactome_pathway_ids(outdir=None):
  """
  Return a dict mapping pathway identifier to pathway title
  Return None if request fails
  """
  rv = None

  # messy logic to use cached version of ReactomePathways.txt
  ofh = None
  line_iter = None
  if outdir is not None:
    fp = os.path.join(outdir, 'ReactomePathways.txt')
    if os.path.exists(fp):
      line_iter = open(fp, 'r')
    else:
      ofh = open(fp, 'w')
  if line_iter is None:
    uri = 'https://reactome.org/download/current/ReactomePathways.txt'
    resp = requests.get(uri)
    if resp.status_code == 200 and len(resp.text) != 0:
      line_iter = resp.iter_lines(decode_unicode=True)

  if line_iter is not None:
    rv = {}
    for line in line_iter:
      line = line.rstrip('\n')
      if ofh is not None:
        ofh.write(line + '\n')
      pathway_id, title, organism = line.split('\t')
      if organism == 'Homo sapiens':
        rv[pathway_id] = title
  else:
    rv = None
  return rv

def get_pathways(pathway_ids, outdir, id_type='reactome', verbose=True):
  """
  Create a file in <outdir> for each pathway in <pathway_ids>
  """
  if id_type != 'reactome' and id_type != 'kegg':
    raise prmf.FactorLibException("Invalid id_type={}; must be one of \"kegg\", \"reactome\".".format(id_type))
  success = []
  fail = []
  last_request = False
  for pathway_id in pathway_ids:
    if verbose:
      sys.stdout.write("Downloading {}...".format(pathway_id))
    pathway_uri = None
    if id_type == 'reactome':
      pathway_uri = reactome_uri(pathway_id)
    else:
      pathway_uri = kegg_uri(pathway_id)

    ofp = os.path.join(args.outdir, pathway_id + ".sif")
    if os.path.exists(ofp):
      if verbose:
        sys.stdout.write("cached\n")
      last_request = False
    else:
      if last_request is True:
        time.sleep(0.5)
      sif_data = pc_get(pathway_uri)
      last_request = True
      if sif_data is not None:
        with open(ofp, 'w') as ofh:
          ofh.write(sif_data)
        success.append(pathway_id)
        sys.stdout.write("success\n")
      else:
        fail.append(pathway_id)
        sys.stdout.write("fail\n")
  return success, fail

def pc_get(uri):
  """
  Return SIF data from uri as a str; return None if request fails
  """
  rv = None
  params= {
    "uri": uri,
    "format": "SIF"
  }
  uri = PC_GET_URI + "?" + urlencode(params)
  resp = requests.get(uri)
  # pathway commons returns 200 even if they arent going to send any data
  if resp.status_code == 200 and len(resp.text) != 0:
    rv = resp.text
  else:
    sys.stderr.write("Request at {} failed; status code: {}\n".format(uri, resp.status_code))
  return rv

def pc_reactome_top_pathways():
  """
  Return a list of pathway identifiers from Reactome from a PathwayCommons TOP_PATHWAYS request
  http://www.pathwaycommons.org/pc2/#top_pathways ; return None if request fails
  """
  rv = []
  uri = 'http://www.pathwaycommons.org/pc2/top_pathways.json'
  params = {
    'q': '*',
    'datasource': 'reactome'
  }
  uri_and_query = '{}?{}'.format(uri, urlencode(params))
  resp = requests.get(uri_and_query)
  if resp.status_code == 200 and len(resp.text) != 0:
    resp_json = resp.json()
    for search_hit in resp_json['searchHit']:
      rv.append(search_hit['uri'])
  else:
    rv = None
  return rv

if __name__ == "__main__":
  parser = argparse.ArgumentParser(description="""
Download Reactome pathways from PathwayCommons in the SIF file format.

While KEGG pathways are nominally present in PathwayCommons, only 1/3 of them are actually available.
As a result, this script does not enable the download of them as is (even though functions to do so are present).
""")
  parser.add_argument("--verbose", "-v", action='store_true', help="If set, print pathway download progress")
  parser.add_argument("--outdir", "-o", help="Directory to save pathways in")
  args = parser.parse_args()

  pathway_to_title = get_reactome_pathway_ids(args.outdir)
  print(len(pathway_to_title))
  ofp = os.path.join(args.outdir, "pathway_titles.tsv")
  with open(ofp, 'w') as fh:
    lines = map(lambda x: '\t'.join(x), pathway_to_title.items())
    fh.write('\n'.join(lines))
  get_pathways(pathway_to_title.keys(), args.outdir, id_type='reactome', verbose=args.verbose)
