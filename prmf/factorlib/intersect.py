"""
Functions related to analysis of intersection of manifold regularizors and iterative updates relating 
to the intersections
"""

import numpy as np
import matplotlib.pyplot as plt

def is_critical_update(x1, x2, c):
  """
  Return True if x2 is on the other side of a critical point c than x1
  """
  dim = x1.shape[0]
  delta = x1 - x2
  slopes = np.ones(dim)
  for i in range(1,dim):
    slopes[i] = delta[i] / delta[0]

  slopes_neg_recip = -1 * np.divide(np.ones(dim), slopes)

  # return vector of length <dim>
  def perp_line(t):
    return c + (t - c[0]) * slopes_neg_recip

  x1_perp = perp_line(x1[0])
  x2_perp = perp_line(x2[0])

  x1_bool = (x1 < x1_perp)
  x2_bool = (x2 < x2_perp)
  x1_and_x2 = np.logical_and(x1_bool, x2_bool)
  rv = np.all(x1_and_x2)

  return rv
