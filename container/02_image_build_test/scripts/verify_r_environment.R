args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: verify_r_environment.R <target-library> <renv.lock>")

target_lib <- args[[1L]]
lockfile <- args[[2L]]
if (!requireNamespace("renv", quietly = TRUE)) {
  stop("renv is not available; set R_LIBS_USER to the bootstrap library")
}

lock <- renv::lockfile_read(lockfile)
locked <- vapply(lock$Packages, function(package) package$Version, character(1))
libraries <- unique(c(target_lib, .Library.site, .Library))
installed <- utils::installed.packages(lib.loc = libraries)
missing <- setdiff(names(locked), rownames(installed))
common <- intersect(names(locked), rownames(installed))
mismatch <- common[locked[common] != installed[common, "Version"]]

cat("R version:", R.version.string, "\n")
cat("Libraries checked:", paste(libraries, collapse = ";"), "\n")
cat("Lockfile package count:", length(locked), "\n")
cat("Matched package versions:", length(common) - length(mismatch), "\n")
cat("Missing packages:", length(missing), "\n")
cat("Version mismatches:", length(mismatch), "\n")
if (length(missing)) cat("Missing:", paste(missing, collapse = ", "), "\n")
if (length(mismatch)) {
  for (package in mismatch) {
    cat("Mismatch:", package, "locked=", locked[[package]],
        "installed=", installed[package, "Version"], "\n")
  }
}
if (length(missing) || length(mismatch)) {
  stop("The installed R libraries do not match all package versions in renv.lock")
}
