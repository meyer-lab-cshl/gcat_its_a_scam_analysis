#!/usr/bin/env Rscript

## -----------------------------------------------------------------------
## parse_benchmark_to_csv.R
##
## Parses a gcatlib benchmark JSON file and writes every per-epoch
## measurement (elapsed time, cpu cycles, instructions, branch
## instructions/misses, page faults, context switches, iterations) to a
## single wide-format CSV -- one row per (benchmark case, epoch).
##
## Quirk this script handles: a benchmark run appends a *new* top-level
## JSON object straight onto the file with no separator, so the file can
## hold several concatenated JSON documents (literally "...}{...").
## This is detected and split automatically -- an ordinary single-document
## file works unchanged.
##
## Usage:
##   Rscript parse_benchmark_to_csv.R <input.json> <output.csv>
## -----------------------------------------------------------------------

suppressPackageStartupMessages({
  library(jsonlite)
  library(dplyr)
  library(purrr)
  library(stringr)
  library(readr)
})

`%||%` <- function(a, b) if (is.null(a)) b else a

# ---- 1. Robust reader: split concatenated JSON documents ---------------

split_json_documents <- function(txt) {
  # Several top-level JSON objects can be glued together with no separator
  # ("...}{..."); mark that boundary with a sentinel not otherwise valid in
  # JSON, then split on it.
  txt <- str_replace_all(txt, "\\}\\s*\\{", "}@@JSON_DOC_SPLIT@@{")
  str_split(txt, fixed("@@JSON_DOC_SPLIT@@"))[[1]]
}

read_bench_json <- function(path) {
  raw  <- read_file(path)
  docs <- split_json_documents(raw)
  parsed <- map(docs, ~ fromJSON(.x, simplifyVector = FALSE))

  # Each document has a top-level $results list; flatten across documents.
  results <- flatten(map(parsed, "results"))

  imap_dfr(results, function(res, i) {
    meas <- res$measurements
    if (is.null(meas) || length(meas) == 0) return(NULL)

    meas_df <- map_dfr(meas, as_tibble)
    meas_df$epoch <- seq_len(nrow(meas_df))

    meas_df %>%
      mutate(
        result_id   = i,
        title       = res$title %||% NA_character_,
        name        = res$name %||% NA_character_,
        unit        = res$unit %||% NA_character_,
        complexityN = res$complexityN %||% NA_real_,
        .before = 1
      )
  })
}

# ---- 2. Load input file and write CSV -----------------------------------

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript parse_benchmark_to_csv.R <input.json> <output.csv>",
       call. = FALSE)
}
input_file <- args[1]
out_file   <- args[2]
stopifnot("Input file not found." = file.exists(input_file))

message("Reading: ", input_file)
bench <- read_bench_json(input_file)
stopifnot("No measurements parsed from input file." = nrow(bench) > 0)

write_csv(bench, out_file)
message("Saved ", nrow(bench), " rows to: ", out_file)
