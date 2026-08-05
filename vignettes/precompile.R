# Run manually after a full `_targets.R` pipeline run:
#   Rscript vignettes/precompile.R
#
# Knits report.Rmd.orig -> report.Rmd, baking in the latest pipeline results
# (tar_read() against the local `_targets` store). The committed report.Rmd
# is then static: R CMD build and pkgdown just typeset it, no `_targets`
# store or network access required.
#
# Runs with the working directory set to vignettes/ -- matches how R CMD
# build/pkgdown knit vignettes, so fig.path output lands next to report.Rmd
# instead of at the repo root.
withr::with_dir("vignettes", knitr::knit(
  input = "report.Rmd.orig",
  output = "report.Rmd"
))
