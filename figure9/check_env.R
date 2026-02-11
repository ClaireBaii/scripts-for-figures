message("R_HOME: ", R.home())
message(".libPaths():")
print(.libPaths())
message("Installed packages:")
installed <- installed.packages()[, "Package"]
if ("mixOmics" %in% installed) {
  message("mixOmics IS installed at: ", find.package("mixOmics"))
} else {
  message("mixOmics is NOT installed.")
}
