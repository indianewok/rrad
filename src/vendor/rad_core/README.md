# Vendored RAD core

The embedded demultiplexing, discovery, spatial, and reformatting engines come
from [RAD](https://github.com/indianewok/rad), under its MIT license. See
[LICENSE](LICENSE). [RAD_VENDOR.json](RAD_VENDOR.json) is the authoritative
record of source and resource commits, adapted files, local patches, and hashes.
Dependency licenses are listed in
[THIRD_PARTY_NOTICES.md](../../../inst/THIRD_PARTY_NOTICES.md) and preserved in
`inst/licenses/`.

## Embedding

The package keeps RAD declarations in `rrad_rad_core` and provides R entry
points instead of a command-line parser. The main adaptations are:

- C++17 filesystem support, gzstream/zlib, and a bounded standard-library
  output queue replace the standalone Boost filesystem/I/O/queue dependencies.
- Native errors reach R as exceptions. Worker errors are collected, interrupts
  are polled on the calling thread, and diagnostics use package-owned capture
  streams rather than process-wide console redirection.
- Reads retain input order. Output files and calibration caches are staged,
  checked, and published with rollback under per-prefix locks. Connections
  remain caller-owned and receive output on the R thread.
- Whitelist roles are explicit and independent of catalog size. Packed
  barcodes support 1–64 bases; wider search windows remain strings.
- HMM calibration follows layout pruning. Position-map caches retain both
  strand models, including parameter-only rows for implicit mirrored states.
- Layout/whitelist discovery and spatial assignment use the same embedded
  resource and output handling as demultiplexing.

The manifest records individual algorithm and compatibility patches. Keep
that detail there rather than duplicating a patch history in this README.
SSW is compiled through `src/rad_core_bridge.cpp`; Edlib has its own translation
unit. Unused plotting and auxiliary upstream headers are excluded.

## Resources

Small layouts and whitelists are installed from `inst/rad/resources/`. The
three large 3M whitelist files are listed in `REMOTE_RESOURCES.tsv` with pinned
URLs, byte counts, and SHA-256 hashes. The untracked upstream
`visium_hd_coordinates.csv.gz` is not bundled.

R verifies downloads and supplies a private, read-only per-call cache overlay.
Native code resolves installed resources or that overlay; it does not search
the working directory, source checkout, or standalone `~/.rad` overrides.
Explicit custom paths remain authoritative. Calibration cache reuse requires both the
layout and position-map files and is enabled only when requested.

Use `rad_download_whitelists()` before an offline run. `rrad.cache_dir` or
`RRAD_CACHE_DIR` sets the cache location; `rrad.offline = TRUE` or
`RRAD_OFFLINE=true` disables downloads. Public resource options are documented
in the R help pages.

## Building

`configure` checks C++17 filesystem support, unsigned 128-bit integers, zlib,
and POSIX threads. It adds `-lstdc++fs` only when needed. Supported 64-bit GCC
and Clang toolchains supply the integer type used by long barcodes.

OpenMP is enabled only after a compile/load probe with R's selected compiler.
`RRAD_OPENMP=auto` is the default; macOS defaults to serial execution.
`RRAD_OPENMP=no` forces serial execution and `RRAD_OPENMP=yes` requires a
working OpenMP toolchain. Results report requested and effective thread
counts, and execution restores the caller's previous OpenMP thread setting.
Boost's supported `BOOST_NO_AUTO_PTR` setting omits obsolete C++17 overloads
without disabling compiler warnings.

## Hashes

`embedded_source_digest` hashes the final patched `include/` tree. Enumerate
its files, sort relative POSIX paths bytewise (`LC_ALL=C`), and build a UTF-8
manifest with one line per file:

```text
<lowercase SHA-256><two spaces><relative path><LF>
```

Hash the complete manifest with SHA-256. Files outside `include/`, including
`RAD_VENDOR.json`, are excluded to avoid self-reference.

`resource_tree_sha256` uses the same procedure relative to
`inst/rad/resources/`. Including `REMOTE_RESOURCES.tsv` commits the hash to
both bundled data and the identities of remote payloads.
