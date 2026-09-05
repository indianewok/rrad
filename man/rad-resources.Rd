\name{rad-resources}
\alias{rad-resources}
\alias{rad_resource_path}
\alias{rad_download_whitelists}
\alias{rad_clear_resource_cache}
\title{Download RAD barcode whitelists on first use}
\usage{
rad_resource_path(key, download = TRUE, cache_dir = NULL, quiet = FALSE)

rad_download_whitelists(kits = NULL, force = FALSE, quiet = FALSE,
  cache_dir = NULL)

rad_clear_resource_cache(confirm = FALSE, cache_dir = NULL)
}
\arguments{
  \item{key}{A built-in whitelist name, RAD kit alias, or whitelist
  filename.}

  \item{download}{Download and save a missing whitelist. If \code{FALSE}, an
  unavailable whitelist returns \code{NA_character_} without accessing the
  network.}

  \item{cache_dir}{Optional cache base directory. The default is
  \code{getOption("rrad.cache_dir")}, then \code{RRAD_CACHE_DIR}, when set,
  otherwise \code{tools::R_user_dir("rrad", "cache")}.}

  \item{quiet}{Suppress download progress messages.}

  \item{kits}{Optional character vector of whitelist names or kit aliases.
  Downloadable names are \code{"10x_3v3"}, \code{"10x_3v4"}, and
  \code{"10x_5v3"}; \code{"10x_3v3.1"} and \code{"10x_3HTv3.1"} are
  accepted aliases for \code{"10x_3v3"}. \code{NULL} selects every
  downloadable whitelist.}

  \item{force}{Re-download and verify whitelists even when valid saved copies
  already exist. A failed refresh leaves the valid saved copy untouched.}

  \item{confirm}{Must be \code{TRUE} to remove downloaded whitelist files.}
}
\value{
\code{rad_resource_path()} returns one normalized file path, or
\code{NA_character_} when \code{download = FALSE} and the whitelist is
unavailable. \code{rad_download_whitelists()} returns a named character vector
of paths.
\code{rad_clear_resource_cache()} invisibly returns whether a cache was
removed.
}
\description{
Most users do not need to call these functions. \code{rad_demux()} and
\code{rad_scan_wl()} automatically download a large built-in barcode
whitelist the first time a kit needs it, then reuse the saved copy on later
runs. The download is written to R's standard per-user cache and checked
against the exact byte count and SHA-256 digest shipped with \code{rrad}.
Package installation, loading, examples, and tests never require network
access.
}
\details{
Set \code{options(rrad.offline = TRUE)} or \code{RRAD_OFFLINE=true} to prohibit
downloads. Already downloaded whitelists remain available. Set
\code{options(rrad.cache_dir = "/shared/path")} or \code{RRAD_CACHE_DIR} to
use a shared or project cache, which can be useful on computing clusters.
}
\examples{
# Locate a saved whitelist without making a network request.
rad_resource_path("10x_3v3", download = FALSE)

\dontrun{
# Normal use is automatic: rad_demux() and rad_scan_wl() download a
# required large built-in whitelist once and reuse it on later calls.
whitelist <- rad_resource_path("10x_3v3")
rad_download_whitelists(c("10x_3v4", "10x_5v3"))

# Optionally download all large whitelists before working offline.
rad_download_whitelists()
}
}
