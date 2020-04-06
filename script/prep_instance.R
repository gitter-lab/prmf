#!/usr/bin/env Rscript
local({
  r <- getOption("repos")
  r["CRAN"] <- "http://cran.r-project.org" 
  options(repos=r)
})
if (!requireNamespace("argparse", quietly = TRUE)) {
  install.packages("argparse")
}
