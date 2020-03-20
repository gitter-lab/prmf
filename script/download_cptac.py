#!/usr/bin/env python
import sys
from urllib.parse import urljoin
import requests
from bs4 import BeautifulSoup
import subprocess as sp

# print to stdout all the protein report summary links from CPTAC
# TODO download the data in the links

PHASE_II_URL = "https://cptc-xfer.uis.georgetown.edu/publicData/Phase_II_Data/"
PHASE_III_URL = "https://cptc-xfer.uis.georgetown.edu/publicData/Phase_III_Data/"

def get_abs_url_links(url):
  sys.stderr.write('Getting links from ' + url + '\n')
  html_doc = requests.get(url).text
  soup = BeautifulSoup(html_doc, 'html.parser')
  abs_url_links = map(lambda x: urljoin(url, x.get('href')), soup.find_all('a'))
  return abs_url_links

def get_protein_summary_for_phase(url):
  study_links = get_abs_url_links(url)
  for cur_url in study_links:
    intra_study_links = get_abs_url_links(cur_url)
    for cur_url_2 in intra_study_links:
      if 'Protein_Report' in cur_url_2:
        protein_report_links = get_abs_url_links(cur_url_2)
        for cur_url_3 in protein_report_links:
          if 'summary.tsv' in cur_url_3:
            print(cur_url_3)

if __name__ == "__main__":
  get_protein_summary_for_phase(PHASE_II_URL)
  get_protein_summary_for_phase(PHASE_III_URL)
