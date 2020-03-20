#!/usr/bin/env Rscript
load_success = TRUE
load_success = load_success && library(BiRewire, quietly=TRUE, logical.return=TRUE)
load_success = load_success && library(igraph, quietly=TRUE, logical.return=TRUE)
load_success = load_success && library(argparse, quietly=TRUE, logical.return=TRUE)
if(!load_success) {
  write("Unable to load dependencies", stderr())
  quit(status=2)
}

# Split file extension off filepath or filename
# TODO copied from PLIER_wrapper.R
splitext = function(text) {
  rv = NA
  pattern = "([^.]*)([.][^.]*)$"
  match_data = regmatches(text, regexec(pattern, text))
  if(length(match_data[[1]]) != 1) {
    rv = list(match_data[[1]][2], match_data[[1]][3])
  }
  rv
}

check_file_and_dir_args = function(args) {
  do_file = FALSE
  do_dir = FALSE
  do_infiles_outdir = FALSE

  if (!is.null(args$infiles) && !is.null(args$outfiles)) {
    do_file = TRUE
    if (length(args$infiles) != length(args$outfiles)) {
      stop('Must provide the same number of infiles and outfiles')
    }
  }
  if (!do_file) {
    if (!is.null(args$indir) && !is.null(args$outdir)) {
      do_dir = TRUE
    }
  }
  if (!do_file && !do_dir) {
    if (!is.null(args$infiles) && !is.null(args$outdir)) {
      do_infiles_outdir = TRUE
    }
  }
  if (!do_file && !do_dir && !do_infiles_outdir) {
    stop('Must provide --infiles and --outfiles or provide --indir and --outdir or provide --infiles and --outdir')
  }

  io_pairs = vector(mode='list')
  if (do_file) {
    io_pairs = Map(c, args$infiles, args$outfiles)
  } else if (do_dir) {
    i = 1
    for (ifn in list.files(path=args$indir)) {
      ifp = file.path(args$indir, ifn)
      ofp = file.path(args$outdir, ifn)
      io_pairs[[i]] = c(ifp, ofp)
      i = i + 1
    }
  } else {
    i = 1
    for (ifp in args$infiles) {
      ofp = file.path(args$outdir, basename(ifp))
      io_pairs[[i]] = c(ifp, ofp)
      i = i + 1
    }
  }
  io_pairs
}

add_file_and_dir_args = function(parser) {
  parser$add_argument('--infiles', nargs='+')
  parser$add_argument('--outfiles', nargs='+')
  parser$add_argument('--indir')
  parser$add_argument('--outdir')
}


main = function() {
  parser = ArgumentParser(description='Randomize edges in a network in a degree-controlled manner')
  add_file_and_dir_args(parser)
  args = parser$parse_args()
  io_pairs = check_file_and_dir_args(args)

  for (i in 1:length(io_pairs)) {
    ifp = io_pairs[[i]][1]
    is_graphml = FALSE
    splitext_rv = splitext(basename(ifp))
    if (!is.na(splitext_rv)) {
      bn = splitext_rv[1]
      ext = splitext_rv[2]
      if (ext == '.graphml') {
        is_graphml = TRUE
      }
    }

    if (is_graphml) {
      ofp = io_pairs[[i]][2]

      G = read_graph(ifp, format='graphml')
      H = birewire.rewire.undirected(G)
      write_graph(H, ofp, format='graphml')
    }
  }

  # also write a list of the outfiles into a single file
  # TODO make this an option
  if (!is.na(args$outdir)) {
    pathways_file = file.path(args$outdir, 'pathways_file.txt')
    for (i in 1:length(io_pairs)) {
      write(io_pairs[[i]][2], append=TRUE, file=pathways_file)
    }
  }
}

main()
