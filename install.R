options(repos = c(CRAN = "https://cloud.r-project.org"))

# Install CRAN packages
install.packages(c("microeco", "file2meco"))

# Optional: fallback to GitHub if CRAN mirror hiccups
if (!requireNamespace("microeco", quietly = TRUE)) {
  if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
  remotes::install_github("ChiLiubio/microeco")
}
if (!requireNamespace("file2meco", quietly = TRUE)) {
  if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
  remotes::install_github("ChiLiubio/file2meco")
}
