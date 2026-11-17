---
name: rrad
description: "Use the rrad R package to inspect, filter, demultiplex, discover barcodes or layouts, and reformat sequencing reads for SIRUP and related workflows. Applies to RRAD analysis tasks and reproducible R scripts, not general package maintenance or standalone C++ RAD development."
---

# RRAD

Use the installed `rrad` R functions. They run the embedded RAD C++ engine;
no separate `rad` executable or service is needed.

## Establish the available API

Run R code in one session, or save a script and run it with `Rscript --vanilla`.
Check the installed package once before choosing a workflow:

```r
library(rrad)
packageVersion("rrad")
find.package("rrad")
rad_core_info()
```

Discovery, explicit spatial export, and HMM options require rrad 1.3.0 or
later. Check `getNamespaceExports("rrad")` and `formals()` when the installed
version differs. If the package is missing or too old, use the user's supplied
checkout or source archive and the package's installation instructions.
Do not assume an unqualified GitHub install contains development-only features.

Read installed help for the operation being used instead of guessing arguments:

```r
help_db <- tools::Rd_db("rrad")
tools::Rd2txt(help_db[["rad-core.Rd"]])  # rad_demux and rad_core_info
```

Other workflow topics are `rad-scan-wl.Rd`, `rad-discover.Rd`,
`rad-spatial-reformat.Rd`, and `rad-reformat.Rd`. Use `names(help_db)` to find
packed-read topics. These manuals are the API reference; this skill does not
duplicate their parameter lists.

## Select the operation

| User's task | R interface |
| --- | --- |
| Inspect or filter a small read file | `rad_read_fastq()`, quality/length columns, `rad_materialize()` |
| Process a large read file in chunks | `rad_stream()`, `rad_stream_next()`, `rad_stream_close()` |
| Compare known barcode catalogs or infer a layout | `rad_discover()` |
| Call observed barcodes from an adapter or capture pattern | `rad_scan_wl()` |
| Demultiplex a specified library | `rad_demux()` |
| Export pRESTO FASTA, change tags, or split by barcode | `rad_reformat()` |
| Assign reads to an experimental spatial map | `rad_spatial_reformat()` |

Use chemistry, adapter sequences, whitelist references, and spatial maps from
the user or their experimental metadata. Catalog matching alone does not
establish chemistry. If those inputs are unresolved, inspect a bounded sample
or report discovery's alternatives before choosing a demultiplexing layout.

`rad_discover()$success` means execution succeeded. Check `complete`,
`layout_complete`, `whitelist_complete`, and match status before using inferred
results. Ties, no hits, and partial layouts are valid outcomes. An explicit
local `catalogs` vector restricts the comparison and avoids default catalog
downloads. Numeric packed catalogs need their declared barcode lengths.

Keep the chosen biological filters and defaults unless the analysis calls for
a change. HMM segmentation is opt-in (`concat_hmm = TRUE`); unresolved parent
reads are not additional demultiplexed molecules. In masked scans, uppercase
`N` captures and lowercase `n` skips. Ordinary scans retain their separate
`after_adapter` behavior. Spatial positions require the actual map or axis
tables and complete `CB:Z:`/`UB:Z:` tags; do not infer positions from whitelist
order or invent missing UMIs. Successful execution can still produce zero
calls: inspect counts and diagnostics before describing a run as recovered data.

## Preserve reads and outputs

- Packed read columns hold native pointers. Keep them in the same R process,
  or use `rad_save()` and `rad_load()` across sessions. Plain `saveRDS()` cannot
  preserve a packed table. `rad_rbind_reads()` combines independent stores;
  direct binding of their handle columns is not safe. Materialize ID keys for
  exact joins between independent stores.
- Filter using the `length` and `mean_q` columns, or the documented quality
  functions, before decoding large sequences. Close streams and connections
  with `on.exit()` inside the function that opens them.
- Choose a fresh output directory or prefix unless replacement was requested.
  Demux and scan can replace their own earlier outputs. Reformat can replace
  its input when `output_fastq` is omitted: pass an explicit new path to retain
  the source. Spatial export requires a new directory.
- Demux, scan, and reformat return `artifacts`; discovery and explicit spatial
  export return `files`. Use those returned paths. A streamed demux has no
  primary file; bytes written to its connection cannot be rolled back.
- Respect offline settings. Large built-in whitelists may need a verified
  download. Use `rad_resource_path(key, download = FALSE)` to check availability
  and `rad_download_whitelists()` when downloads are part of the task.
- Use the requested thread budget and inspect the effective count. A serial
  build uses one worker even when more were requested.

Deliver the requested data products and a short report of input/output counts,
filter settings, discovery uncertainty or unresolved reads, and output paths.
Retain the returned run configuration, core provenance, and logs with the
analysis. Validate the produced records or summaries relevant to the task;
package-wide CRAN and regression runs belong to package maintenance, not every
analysis.
