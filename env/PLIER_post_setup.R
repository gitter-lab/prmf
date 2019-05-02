#!/usr/bin/env Rscript
library(devtools)
repos = getOption("repos")
repos["CRAN"] = "http://cran.us.r-project.org"
options(repos=repos)
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("qvalue")
# ref the latest commit before PLIERconditional (which breaks install)
install_github("wgmao/PLIER", ref="afb4ccbf761418535c3e47ec6baeb6bcd8ec716a")
