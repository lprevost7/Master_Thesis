#!/usr/bin/env Rscript

# Usage:
#   Rscript hypervolume_from_csv.R graph.csv

suppressPackageStartupMessages({
  library(moocore)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript hypervolume_from_csv.R path/to/graph.csv")
}

csv_file <- args[1]

df <- read.csv(csv_file)

points <- as.matrix(df[, c("X", "Y")])

points <- unique(points)

points_nd <- moocore::filter_dominated(points)

ref <- apply(points_nd, 2, max) * 1.01

hv <- moocore::hypervolume(points_nd,reference = ref,maximise = c(FALSE, FALSE))

cat(hv, "\n")
