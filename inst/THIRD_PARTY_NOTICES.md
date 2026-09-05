# Third-party notices

rrad is licensed under the GNU General Public License, version 3 or later.
That package license does not replace the licenses and notices below.

This inventory covers the vendored native include closure used by rrad. In
particular, `src/rad_core_bridge.cpp` includes RAD's `rad/rad_headers.h`, which
transitively includes the libraries listed here. Dedicated package translation
units compile the vendored Edlib implementation, and the RAD bridge textually
includes the SSW C and C++ implementation files. `sse2neon.h` is used by SSW
on ARM/NEON builds. The corresponding complete license texts are installed in
`licenses/`.

## RAD embedded core

- Component: RAD core headers and rrad's embedded-core adaptations
- Vendored paths: `src/vendor/rad_core/include/rad/*`
- License: MIT
- Copyright: Copyright (c) 2025 cmv
- License text: `licenses/RAD-MIT.txt`

## RAD runtime data resources

The read layouts and smaller whitelist tables installed beneath
`rad/resources/` come from the version-pinned RAD source tree identified in
`rad/RAD_RESOURCES.json`. Three larger 3M whitelist tables are not shipped in
the rrad package archive. Their immutable upstream commit URLs, exact byte
counts, SHA-256 digests, and kit aliases are recorded in
`rad/resources/REMOTE_RESOURCES.tsv`; rrad retrieves a required table at
runtime and uses it only after integrity verification.

The URL and digest catalog records provenance and integrity. It does not
replace or modify any terms that apply to the upstream data.

## CSV Parser and its embedded components

- Component: CSV Parser for C++, version 2.3.0
- Vendored path: `src/vendor/rad_core/include/csv-parser/single_include/csv.hpp`
- License: MIT
- Copyright: Copyright (c) 2017-2024 Vincent La
- License text: `licenses/csv-parser-MIT.txt`

The single-header CSV Parser copy embeds the following separately attributed
code:

- `mio` memory-mapping code: MIT; Copyright 2017
  https://github.com/mandreyel. See `licenses/mio-MIT.txt`.
- `string-view-lite`, version 1.1.0: Boost Software License 1.0; Copyright
  2017-2019 Martin Moene. See `licenses/Boost-1.0.txt`.
- Hedley, version 9: CC0 1.0 Universal; created by Evan Nemerson. See
  `licenses/CC0-1.0.txt`.
- JSON escaping functions modified from JSON for Modern C++: MIT; Copyright
  2013-2015 Niels Lohmann. See `licenses/nlohmann-json-MIT.txt`.

## Parallel Hashmap

- Component: Parallel Hashmap, version 2.0.0, including modified work from
  Abseil
- Vendored paths used by the include closure:
  `src/vendor/rad_core/include/parallel_hashmap/phmap.h`, `phmap_base.h`,
  `phmap_bits.h`, `phmap_config.h`, `phmap_fwd_decl.h`, and `phmap_utils.h`
- License: Apache License 2.0
- Copyright: Copyright (c) 2019 Gregory Popovitch
- Copyright: Copyright 2018 The Abseil Authors
- License text: `licenses/Apache-2.0.txt`

`phmap_utils.h` also contains hash-combine code under the Boost Software
License 1.0, Copyright 2005-2014 Daniel James. See
`licenses/Boost-1.0.txt`.

## Edlib

- Component: Edlib public C API header and implementation
- Vendored paths: `src/vendor/rad_core/include/edlib/include/edlib.h` and
  `src/vendor/rad_core/include/edlib/src/edlib.cpp`
- License: MIT
- Copyright: Copyright (c) 2014 Martin Šošić
- License text: `licenses/Edlib-MIT.txt`

## gzstream

- Component: gzstream C++ iostream wrappers for zlib
- Vendored path: `src/vendor/rad_core/include/gzstream/gzstream.h`
- License: GNU Lesser General Public License 2.1 or, at the recipient's option,
  any later version
- Copyright: Copyright (C) 2001 Deepak Bandyopadhyay, Lutz Kettner
- License text: `licenses/LGPL-2.1.txt`

## Complete Striped Smith-Waterman Library (SSW)

- Component: SSW C library and C++ wrapper
- Vendored paths: `src/vendor/rad_core/include/ssw/ssw.c`, `ssw.h`,
  `ssw_cpp.cpp`, and `ssw_cpp.h`
- License: MIT, with the BSD-2-Clause material described below
- Copyright: Copyright 2010 Boston College
- Copyright: Copyright (c) 2012-2015 Boston College
- Copyright: Copyright (c) 2012-2025 Mengyao Zhao and other SSW developers
- License text: `licenses/SSW-MIT.txt`

The SSW implementation contains work derived from Michael Farrar's alignment
code under the BSD 2-Clause License, Copyright 2006 Michael Farrar. See
`licenses/SSW-BSD-2-Clause.txt`. Its source also identifies the lazy-F loop as
derived from MIT-licensed SWPS3 work, Copyright (c) 2007-2008 ETH Zürich,
Institute of Computational Science. See `licenses/SWPS3-MIT.txt`.

On ARM/NEON builds, SSW includes `ssw/sse2neon.h`, an MIT-licensed translation
layer. The vendored file credits John W. Ratcliff, Brandon Rowlett, Ken Fast,
Eric van Beurden, Alexander Potylitsin, Hasindu Gamaarachchi, Jim Huang, Mark
Cheng, Malcolm James MacLeod, Devin Hussey, Sebastian Pop, Developer Ecosystem
Engineering, Danila Kutenin, François Turban, Pei-Hsuan Hung, Yang-Hao Yuan,
Syoyo Fujita, and Brecht Van Lommel. See `licenses/sse2neon-MIT.txt`.

## kseq

- Component: kseq FASTA/FASTQ parser
- Vendored path: `src/vendor/rad_core/include/kseq/kseq.h`
- License: MIT
- Copyright: Copyright (c) 2008, 2009, 2011 Attractive Chaos
  <attractor@live.co.uk>
- License text: `licenses/kseq-MIT.txt`

## fkYAML

- Component: fkYAML, version 0.4.2
- Vendored path: `src/vendor/rad_core/include/fkYAML/node.hpp`
- License: MIT
- Copyright: Copyright 2023-2025 Kensuke Fukutani
  <fktn.dev@gmail.com>
- License text: `licenses/fkYAML-MIT.txt`

## Audit boundary

This notice inventory is limited to code vendored beneath
`src/vendor/rad_core/include`. System libraries and headers (for example zlib
and the platform C/C++ runtime) and dependencies supplied by other R packages
(for example Boost headers supplied by BH) are not redistributed as part of
this vendored tree and are outside this inventory.
