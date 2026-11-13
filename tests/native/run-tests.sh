#!/usr/bin/env bash
# Standalone native regression suite; no package source files are compiled in place.
set -euo pipefail
project_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
build_dir="$(mktemp -d "${TMPDIR:-/tmp}/rrad-native-tests.XXXXXX")"
trap 'rm -rf "$build_dir"' EXIT
read -r -a compiler <<< "$(R CMD config CXX17)"
read -r -a standard <<< "$(R CMD config CXX17STD)"
bh_include="$(Rscript -e 'cat(system.file("include", package="BH"))')"
if [[ -z "$bh_include" ]]; then
  echo 'Install the BH R package before running native tests.' >&2
  exit 1
fi
vendor_dir="$project_dir/src/vendor/rad_core/include"
flags=("${standard[@]}" -O1 -pthread -I"$bh_include" -I"$vendor_dir"
       -I"$vendor_dir/edlib" -I"$vendor_dir/edlib/include" -I"$vendor_dir/ssw"
       "-DRAD_SOURCE_DIR=\"$project_dir/inst/rad\"")
# Set RRAD_NATIVE_FLAGS for optional sanitizers/OpenMP; keep default tests serial.
if [[ -n "${RRAD_NATIVE_FLAGS:-}" ]]; then
  read -r -a extra_flags <<< "$RRAD_NATIVE_FLAGS"
  flags+=("${extra_flags[@]}")
fi
"${compiler[@]}" "${flags[@]}" -x c++ -c "$vendor_dir/ssw/ssw.c" -o "$build_dir/ssw.o"
"${compiler[@]}" "${flags[@]}" -c "$vendor_dir/ssw/ssw_cpp.cpp" -o "$build_dir/ssw_cpp.o"
"${compiler[@]}" "${flags[@]}" -c "$vendor_dir/edlib/src/edlib.cpp" -o "$build_dir/edlib.o"
if [[ $# -eq 0 ]]; then set -- "$project_dir"/tests/native/*_test.cpp; fi
for source in "$@"; do
  name="$(basename "$source" .cpp)"
  echo "Building $name"
  "${compiler[@]}" "${flags[@]}" "$source" "$build_dir/ssw.o" "$build_dir/ssw_cpp.o" \
    "$build_dir/edlib.o" -lz -o "$build_dir/$name"
  "$build_dir/$name"
done
