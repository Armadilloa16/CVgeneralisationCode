

fpath = "../CVgeneralisations/code/data/sim_results"
files = list.files(fpath)
for (f in files) {
  cat(f, '\n')
  tmp = read.csv(file.path(fpath, f))
  tmp$dataset = NULL
  tmp$bin.shift = NULL
  if ('method' %in% names(tmp)) {
    tmp$method = NULL
  }
  write.csv(tmp, file.path(fpath, f), row.names = FALSE)
}