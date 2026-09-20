args <- commandArgs(trailingOnly = FALSE)
script <- sub("^--file=", "", args[startsWith(args, "--file=")][[1]])
root <- dirname(normalizePath(script, winslash = "/", mustWork = TRUE))
expected <- read.csv(file.path(root, "R_packages.csv"))
expected$installed <- vapply(expected$package, function(p) {
  if (!requireNamespace(p, quietly = TRUE)) return("missing")
  as.character(utils::packageVersion(p))
}, character(1))
print(expected, row.names = FALSE)
if (as.character(getRversion()) != "4.5.3" || any(expected$version != expected$installed)) {
  stop("The installed R environment differs from the recorded analysis environment.")
}
message("R and all seven direct dependency versions match.")
