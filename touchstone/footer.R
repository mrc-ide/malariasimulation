# You can modify the PR comment footer here. You can use github markdown e.g.
# emojis like :tada:.
# This file will be parsed and evaluate within the context of
# `benchmark_analyze` and should return the comment text as the last value.
# See `?touchstone::pr_comment`

documentation <- "https://lorenzwalthert.github.io/touchstone/articles/inference.html"

# This is exported by the workflow itself
workflow <- Sys.getenv("WORKFLOW_URL")

glue::glue(
  "\n\nFurther explanation regarding interpretation and",
  " methodology can be found in the [documentation]({documentation}).",
  "\nPlots and raw data are available as artifacts of",
  " [the workflow run]({workflow})."
)
