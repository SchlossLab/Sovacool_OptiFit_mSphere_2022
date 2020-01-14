get_contigsfile <- function(){
  p <- "data/mouse/raw/"
  p <- gsub("\\/$", "", p)
  f <- list.files(path=p, pattern="*.fastq.gz")
  f <- f[!grepl("Mock", f)] #let's ignore the mock community data

  r1 <- f[grep("_R1_", f)]
  r1.group <- sub("_S.*", "", r1)

  r2 <- f[grep("_R2_", f)]
  r2.group <- sub("_S.*", "", r2)

  stopifnot(r1.group == r2.group)

  files.file <- paste0(gsub(".*\\/(.*)", "\\1", p), ".files")
  p.files.file <- paste0(p, "/", files.file)
  output_table <- cbind(r1.group, r1, r2)
  write.table(file="data/mouse/mouse.files", output_table, sep="\t", quote=F, row.names=F, col.names=F)
}

get_contigsfile()
