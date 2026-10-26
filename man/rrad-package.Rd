\name{rrad-package}
\alias{rrad-package}
\alias{rrad}
\docType{package}
\title{\packageTitle{rrad}}
\description{\packageDescription{rrad}}
\details{
\pkg{rrad} is the RAD companion package for SIRUP. The embedded RAD engines
run in the R process on Linux and macOS without a separate command-line
installation. Use \code{\link{rad_core_info}()} to inspect native provenance
and build capabilities, \code{\link{rad_demux}()} to demultiplex,
\code{\link{rad_scan_wl}()} to discover a barcode whitelist, and
\code{\link{rad_reformat}()} to transform tagged reads.

\code{\link{rad_pack_reads}()} provides compact read tables;
\code{\link{rad_stream}()} supports chunked FASTQ processing.
When installed with vignettes, \code{vignette("rrad-workflow", package = "rrad")}
opens a complete offline example.

Small RAD layouts and whitelist tables ship with the package. Larger built-in
3M barcode whitelists are downloaded only when a matching kit is first used,
verified against packaged byte counts and SHA-256 digests, and reused on later
runs. Use \code{\link{rad_download_whitelists}()} to download them
ahead of offline or batch work and \code{\link{rad_resource_path}()} to find
a saved whitelist file.
}
\section{Package Content}{\packageIndices{rrad}}
\author{\packageAuthor{rrad}}
\section{Maintainer}{\packageMaintainer{rrad}}
\keyword{package}
