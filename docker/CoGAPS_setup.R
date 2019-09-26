#!/usr/bin/env Rscript
requireAndInstall = function(package) {
  if (!requireNamespace(package, quietly=TRUE)) {
    install.packages(package)
  }
}
requireAndInstall("remotes")
requireAndInstall("argparse")
requireAndInstall("BiocManager")
BiocManager::install("FertigLab/CoGAPS")
