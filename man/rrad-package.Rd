\name{rrad-package}
\alias{rrad-package}
\alias{rrad}
\docType{package}
\title{\packageTitle{rrad}}
\description{\packageDescription{rrad}}
\details{
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
