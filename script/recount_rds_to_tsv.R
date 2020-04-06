#!/usr/bin/env Rscript
requireNamespace('argparse')
main = function() {
  parser = argparse::ArgumentParser(description=paste0("",
"Assumes <--infile> is an RDS file which can be the target of write.table",
"",
"TODO",
"----",
"- Type checking",
""))
  parser$add_argument('--infile', '-i', help="rds file recount_rpkm.RDS from zipped archive at https://ndownloader.figshare.com/files/10881866")
  parser$add_argument('--outdir', '-o', help="Output directory")
  args = parser$parse_args()

  rpkm.df = readRDS(args$infile)

  # remove ENSG version numbers
  # https://github.com/greenelab/rheum-plier-data/blob/978c37938383ff7adcadacfcbc35931ce5e62b17/recount2/2-prep_recount_for_plier.R#L34
  genes = unlist(lapply(strsplit(rpkm.df$ENSG, "[.]"), `[[`, 1))
  samples = colnames(rpkm.df)[-1]
  gene_col_name = colnames(rpkm.df)[1]

  # transpose data to sample x gene
  # https://stackoverflow.com/questions/6778908/transpose-a-data-frame
  rpkm.df = t(rpkm.df[,-1])
  colnames(rpkm.df) = genes

  write.table(rpkm.df, file.path(args$outdir, sub('.rds', '.tsv', basename(args$infile), ignore.case=TRUE)), sep='\t', col.names=TRUE, row.names=TRUE)
  #write.table(rpkm.df[1:5,1:5], file.path(args$outdir, sub('.rds', '.slice.tsv', basename(args$infile), ignore.case=TRUE)), sep='\t', col.names=TRUE, row.names=TRUE)
}

main()
