#!/usr/bin/env python

from distutils.core import setup

setup(name="prmf",
  version="0.1",
  description="Library and scripts for Pathway-Regularized Matrix Factorization",
  author="Aaron Baker",
  author_email="abaker@cs.wisc.edu",
  packages=['prmf'],
  scripts=[
    'script/prmf_runner.py'
  ]
)
