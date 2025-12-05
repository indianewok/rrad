#!/usr/bin/env Rscript

cat("Welcome to rad! Checking which packages are necessary...\n")
libraries<-c("docopt", "data.table", "magrittr", "stringr", "stringdist","bit64","mclust","pracma", "ggplot2","ggpubr")
github_libraries<-c("rad" = "indianewok/rad@dev")
install_and_load<-function(lib, repo = NULL) {
  if (!require(lib, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)) {
    if (is.null(repo)) {
      cat(paste0("Installing CRAN package: ", lib, "...\n"))
      install.packages(lib, repos = "http://cran.us.r-project.org", dependencies = TRUE, quiet = TRUE)
    } else {
      if (!require(remotes, quietly = TRUE, warn.conflicts = FALSE)) {
        cat("Installing CRAN package: remotes...\n")
        install.packages("remotes", repos = "http://cran.us.r-project.org", quiet = TRUE)
      }
      cat(paste0("Installing GitHub package: ", lib, " from ", repo, "...\n"))
      remotes::install_github(repo, quiet = TRUE)
    }
    cat(paste0("Loading package: ", lib, "...\n"))
    invisible(library(lib, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE))
  } else {
    cat(paste0("Package ", lib, " is already installed and loaded.\n"))
  }
}
cat("Checking CRAN packages...\n")
invisible(lapply(libraries, install_and_load))
cat("Checking GitHub packages...\n")
invisible(lapply(names(github_libraries), function(lib){
  install_and_load(lib, github_libraries[lib])
}))
cat("All necessary R packages are installed and loaded.\n")


doc <- "
Usage: rad_run.R [--fastq_file_or_directory_path=<path>] [--read_layout_path=<path>] [--output_directory_path=<path>] [--compress] [--generate_sigstring_diagnostics] [--sigstring_diagnostic_verbose] [--misalignment_threshold_path=<path>] [--tabulated_sigstring_count=<count>] [--chunk_size=<size>] [--nthreads=<threads>]

Options:
  -h --help                                Show this screen.
  -i --fastq_file_or_directory_path=<path> Path to your input file(s).
  -r --read_layout_path=<path>             Path to your read_layout csv.
  -o --output_directory_path=<path>        Path to your output directory, if it doesn't exist I'll make one.
  -c --compress                            Compress the output (default: TRUE, compresses to .gz).
  -d --generate_sigstring_diagnostics      Generate sigstring diagnostics, a summary table specifically. (default: FALSE).
  -v --sigstring_diagnostic_verbose        Generate verbose sigstring diagnostics for all sigstrings. Not recommended for super-large datasets 'cuz it takes too long (default: FALSE).
  -m --misalignment_threshold_path=<path>  Path to the misalignment threshold file. Use if you want to play with the thresholds to see if more flexible parameters will help read recovery (default: NULL).
  -s --tabulated_sigstring_count=<count>   Number of sigstrings to tabulate if generating diagnostics is TRUE, default is all [default: -1].
  -z --chunk_size=<size>                   Size of the chunks for processing [default: 50000].
  -t --nthreads=<threads>                  Number of threads to use [default: 1].
  "
  
  arguments<-docopt(doc)
  
  # Convert arguments to the appropriate types
  arguments$`--tabulated_sigstring_count` <- as.integer(arguments$`--tabulated_sigstring_count`)
  arguments$`--chunk_size`<-as.integer(arguments$`--chunk_size`)
  arguments$`--nthreads`<-as.integer(arguments$`--nthreads`)
  arguments$`--compress`<-!is.null(arguments$`--compress`)
  arguments$`--generate_sigstring_diagnostics`<-!is.null(arguments$`--generate_sigstring_diagnostics`)
  arguments$`--sigstring_diagnostic_verbose`<-!is.null(arguments$`--sigstring_diagnostic_verbose`)
  
  # Ensure the function is defined
  if (!exists("rad_run", where = "package:rad", mode = "function")) {
    stop("rad_run function is not available in the rad package.")
  }
  
  # Print the values to ensure they are correctly parsed
  cat("Calling rad_run with the following parameters:\n")
  cat("fastq_file_or_directory_path:", arguments$`--fastq_file_or_directory_path`, "\n")
  cat("read_layout_path:", arguments$`--read_layout_path`, "\n")
  cat("output_directory_path:", arguments$`--output_directory_path`, "\n")
  cat("compress:", arguments$`--compress`, "\n")
  cat("generate_sigstring_diagnostics:", arguments$`--generate_sigstring_diagnostics`, "\n")
  cat("sigstring_diagnostic_verbose:", arguments$`--sigstring_diagnostic_verbose`, "\n")
  cat("misalignment_threshold_path:", arguments$`--misalignment_threshold_path`, "\n")
  cat("tabulated_sigstring_count:", arguments$`--tabulated_sigstring_count`, "\n")
  cat("chunk_size:", arguments$`--chunk_size`, "\n")
  cat("nthreads:", arguments$`--nthreads`, "\n")
  
  rad_run(fastq_file_or_directory_path = arguments$fastq_file_or_directory_path,
    read_layout_path = arguments$`--read_layout_path`,
    output_directory_path = arguments$`--output_directory_path`,
    compress = arguments$`--compress`,
    generate_sigstring_diagnostics = arguments$`--generate_sigstring_diagnostics`,
    sigstring_diagnostic_verbose = arguments$`--sigstring_diagnostic_verbose`,
    misalignment_threshold_path = arguments$`--misalignment_threshold_path`,
    tabulated_sigstring_count = arguments$`--tabulated_sigstring_count`,
    chunk_size = arguments$`--chunk_size`,
    nthreads = arguments$`--nthreads`
  )
