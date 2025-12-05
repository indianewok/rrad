#!/usr/bin/env Rscript

library(docopt)
library(rad)
doc <- "
Usage:
  rad_cli_run.R --read_layout_form=<rlf> --read_dir=<rd> --output_dir=<od> --whitelist_path=<wp> [--nthreads=<nt>] [--SETUP_DRY_RUN] [--FULL_RUN] [--VERBOSE] [--KINDA_VERBOSE] [--ORIGINAL_CHUNK_SIZE=<ocs>]

Options:
  --read_layout_form=<rlf>       Read layout form file path.
  --read_dir=<rd>                Directory containing read files.
  --output_dir=<od>              Directory where output will be stored.
  --whitelist_path=<wp>          Path to whitelist file.
  --nthreads=<nt>                Number of threads to use [default: 1].
  --SETUP_DRY_RUN                Perform a setup dry run.
  --FULL_RUN                     Perform a full run.
  --VERBOSE                      Enable verbose output.
  --KINDA_VERBOSE                Enable kind-of-verbose output.
  --ORIGINAL_CHUNK_SIZE=<ocs>    Original chunk size [default: 500000].
"

# Parse command-line arguments
args<-docopt(doc)

args$SETUP_DRY_RUN<-!is.null(args[['--SETUP_DRY_RUN']])
args$FULL_RUN<-!is.null(args[['--FULL_RUN']])
args$VERBOSE<-!is.null(args[['--VERBOSE']])
args$KINDA_VERBOSE<-!is.null(args[['--KINDA_VERBOSE']])
args$nthreads<-as.integer(args[['--nthreads']])
args$ORIGINAL_CHUNK_SIZE<-as.numeric(args[['--ORIGINAL_CHUNK_SIZE']])

rad::rad_run(
  read_layout_form = args[['--read_layout_form']],
  read_dir = args[['--read_dir']],
  output_dir = args[['--output_dir']],
  whitelist_path = args[['--whitelist_path']],
  nthreads = args$nthreads,
  SETUP_DRY_RUN = args$SETUP_DRY_RUN,
  FULL_RUN = args$FULL_RUN,
  VERBOSE = args$VERBOSE,
  KINDA_VERBOSE = args$KINDA_VERBOSE,
  ORIGINAL_CHUNK_SIZE = args$ORIGINAL_CHUNK_SIZE
)