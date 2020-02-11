---
title: 'Pathway-Regularized Matrix Factorization'
tags:
  - Python
  - matrix factorization
  - graph
  - network
authors:
  - name: Aaron Baker
    orcid: 0000-0002-2815-9932
    affiliation: '1, 2' # (Multiple affiliations must be quoted)
  - name: Anthony Gitter
    orcid: 0000-0002-5324-9833
    affiliation: '1, 2, 3'
affiliations:
 - name: Department of Computer Sciences, University of Wisconsin-Madison
   index: 1
 - name: Morgridge Institute for Research
   index: 2
 - name: Department of Biostatistics and Medical Informatics, University of Wisconsin-Madison
   index: 3
date: 11 February 2020
bibliography: paper.bib
---

# Summary

Pathway-Regularized Matrix Factorization (PRMF) is an extension of non-negative matrix factorization for high-throughput biological data.
It uses graph structures from biological pathways to constrain the factorization.

![Pathway-Regularized Matrix Factorization overview.](overview.png)

Related work:
- Review [@stein-obrien_enter_2018]
- pyNBS [@huang_pynbs_2018]
- PLIER [@mao_plier_2019]

PRMF name disambiguation:
- Probabilistic relational matrix factorization [@liu_prmf_2016]
- Probabilistic robust matrix factorization [@hutchison_prmf_2012]

# Mathematics

Single dollars ($) are required for inline mathematics e.g. $f(x) = e^{\pi/x}$

Double dollars make self-standing equations:

$$\Theta(x) = \left\{\begin{array}{l}
0\textrm{ if } x < 0\cr
1\textrm{ else}
\end{array}\right.$$


# Citations

Citations to entries in paper.bib should be in
[rMarkdown](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html)
format.

For a quick reference, the following citation commands can be used:
- `@author:2001`  ->  "Author et al. (2001)"
- `[@author:2001]` -> "(Author et al., 2001)"
- `[@author1:2001; @author2:2001]` -> "(Author1 et al., 2001; Author2 et al., 2002)"

# Acknowledgements

**Add funding and individual acknowledgements**

# References
