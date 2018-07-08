#!/usr/bin/env python
import argparse
import numpy as np
from sklearn.metrics import precision_recall_curve
import random
import matplotlib.pyplot as plt

if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument("outfile")
  args = parser.parse_args()

  n = 1000
  n_pos = 100
  true_ind = random.sample(range(n), n_pos)
  y_true = np.zeros((1000,))
  for ind in true_ind:
    y_true[ind] = 1
  y_score = np.random.rand(1000)

  precision, recall, thresholds = precision_recall_curve(y_true, y_score)

  plt.clf()
  plt.step(recall, precision, color='red', where='post')
  plt.fill_between(recall, precision, step='post', alpha=0.2, color='r')
  plt.xlabel('Recall', fontsize='x-large')
  plt.ylabel('Precision', fontsize='x-large')
  plt.ylim([0.0, 1.05])
  plt.xlim([0.0, 1.0])
  plt.title('Precision-Recall')
  plt.legend()
  plt.savefig(args.outfile, bbox_inches='tight')
