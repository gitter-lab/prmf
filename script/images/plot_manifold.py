#!/usr/bin/env python
import sys, argparse
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import itertools as it
import math

def plot_gallery_manifolds(title, images, manifolds, n_col=3, n_row=2, edge_limit=math.inf):
  image_shape = (64, 64) # TODO
  plt.figure(figsize=(2. * n_col, 2.26 * n_row))
  plt.suptitle(title, size=16)

  def int_to_index(x):
    row = x // image_shape[1] 
    col = x % image_shape[1]
    return row, col

  def plot_manifold(plt, manifold):
    # pairs of x,y coordinates for each line segment
    segs = []
    k = 0
    for i, j in manifold.edges_iter():
      i = int(i)
      j = int(j)
      # convert pixel-pixel coordinates
      row_i, col_i = int_to_index(i)
      row_j, col_j = int_to_index(j)
      k += 1
      if k >= edge_limit:
        break

      # define line segment between coordinates
      segs.append(((col_i, row_i), (col_j, row_j))) # transpose for use with imshow
    coll = LineCollection(segs)
    ax = plt.gca()
    ax.add_collection(coll)

  for i, comp in enumerate(images):
    manifold = manifolds[i]
    plt.subplot(n_row, n_col, i + 1)
    vmax = max(comp.max(), -comp.min())
    plt.imshow(comp.reshape(image_shape), cmap=plt.cm.gray,
           interpolation='nearest',
           vmin=-vmax, vmax=vmax)
    plot_manifold(plt, manifold)
    plt.xticks(())
    plt.yticks(())
  plt.subplots_adjust(0.01, 0.05, 0.99, 0.93, 0.04, 0.)

def main():
  parser = argparse.ArgumentParser(description="""
Plot pixel-pixel distances over images
""")
  parser.add_argument("--components", "-c", help="comma delimited matrix")
  parser.add_argument("--manifolds", "-l", help="one or more files of weighted edge lists", nargs="+")
  parser.add_argument("--edge-limit", default=math.inf, type=int)
  parser.add_argument("--outfile", "-o", help=".png to write")
  args = parser.parse_args()

  comps = np.genfromtxt(args.components, delimiter=",")
  manifolds = list(map(nx.read_graphml, args.manifolds))

  plot_gallery_manifolds("Components", comps, manifolds, edge_limit=args.edge_limit)
  plt.savefig(args.outfile)

if __name__ == "__main__":
  main()
