#!/usr/bin/env python
import argparse
import networkx as nx
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

# https://www.infobyip.com/detectmonitordpi.php
DPI = 96

def main():
  parser = argparse.ArgumentParser()
  parser.add_argument("out", type=str, help="png")
  args = parser.parse_args()

  fig = plt.figure(figsize=(900/DPI, 900/DPI), dpi=DPI)

  ax = fig.gca()
  ax.set_xticks(list(map(lambda x: x / 10.0, range(10))))
  ax.grid(which='both', markevery=0.05)
  #ax.axis('off')
  ax.set_xlim([0.0,1.0])
  #ax.xaxis.set_major_locator(mpl.ticker.NullLocator())
  #ax.yaxis.set_major_locator(mpl.ticker.NullLocator())
  fig.savefig(args.out)
  #fig.savefig(args.out, pad_inches=0, bbox_inches='tight', transparent=True)

if __name__ == "__main__":
  main()
