#!/usr/bin/env python
import sys, argparse
import matplotlib.pyplot as plt
import prmf
import numpy as np

def main():
  parser = argparse.ArgumentParser(description="""Plot Olivetti faces decomposition""")
  parser.add_argument("--delimiter", "-d", default="\t", type=str)
  parser.add_argument("--infile", "-i", help="tsv file of components to visualize")
  parser.add_argument("--outfile", "-o", help="png to write")
  args = parser.parse_args()

  arr = np.genfromtxt(args.infile, delimiter=args.delimiter)
  n_row, n_col = arr.shape
  n_pixel = prmf.OLIVETTI_SHAPE[0] * prmf.OLIVETTI_SHAPE[1]
  images = []
  if(n_row == n_pixel):
    for i in range(n_col):
      images.append(arr[:,i].reshape(prmf.OLIVETTI_SHAPE))
  elif(n_col == n_pixel):
    for i in range(n_row):
      images.append(arr[i,:].reshape(prmf.OLIVETTI_SHAPE))
  else:
    raise RuntimeError("infile has invalid shape: neither {} nor {} matches the number of pixels in the Olivetti faces dataset {}".format(n_row, n_col, n_pixel))
  prmf.plot_olivetti_components(plt, images, title='Components')
  plt.savefig(args.outfile)

if __name__ == "__main__":
  main()
