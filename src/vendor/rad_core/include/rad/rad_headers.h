#pragma once

// ————————————————————————————————————————————————————————————————
// C++ Standard Library in alphabetical order
// ————————————————————————————————————————————————————————————————
#include <algorithm>
#include <array>
#include <cctype>
#include <cerrno>
#include <climits>
#include <chrono>
#include <cmath>
#include <condition_variable>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <exception>
#include <filesystem>
#include <fstream>
#include <functional>
#include <getopt.h>
#include <iomanip>
#include <iostream>
#ifdef __APPLE__
#include <mach/mach.h>
#endif
#include <map>
#include <limits>
#include <memory>
#include <mutex>
#include <numeric>
#include <optional>
#include <queue>
#include <random>
#include <regex>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <system_error>
#include <thread>
#include <tuple>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <sys/resource.h>
#if defined(__linux__)
#include <unistd.h>
#endif

// ————————————————————————————————————————————————————————————————
// Parallelism
// ————————————————————————————————————————————————————————————————
#include <atomic>
#include "openmp_compat.hpp"

// ————————————————————————————————————————————————————————————————
// Boost
// ————————————————————————————————————————————————————————————————
#include <boost/bimap.hpp>
#include <boost/optional.hpp>
#include <boost/functional/hash.hpp>
#include <boost/multi_index/member.hpp>
#include <boost/bimap/unordered_set_of.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index/sequenced_index.hpp>
#include <boost/math/distributions/poisson.hpp>
#include <boost/math/special_functions/fpclassify.hpp>
#include <boost/multi_index/composite_key.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/multi_index_container.hpp>

// ————————————————————————————————————————————————————————————————
// External libraries
// ————————————————————————————————————————————————————————————————
#include "csv-parser/single_include/csv.hpp"
#include "parallel_hashmap/phmap.h"
#include "edlib/include/edlib.h"
#include "gzstream/gzstream.h"
#include "ssw/ssw_cpp.h"
#include "kseq/kseq.h"
#include "zlib.h"

// Keep upstream's normal terminal streams for the CLI, while allowing an
// embedding host to provide a package-safe captured stream.
#include "embedded_console.hpp"

// ————————————————————————————————————————————————————————————————
//  rad headers
// ————————————————————————————————————————————————————————————————
// Keep the embedded engine in a package-owned namespace so upstream global
// declarations cannot collide with other native code loaded into R.
namespace rrad_rad_core {

#include "misc_utils.hpp"
#include "io_streaming.hpp"
#include "barcode_correction.hpp"
#include "adapter_thresholds.hpp"
#include "read_layout.hpp"
#include "sigstring.hpp"
#include "whitelist_generator.hpp"
#include "reformat.hpp"

}  // namespace rrad_rad_core

// seqspec carries a third-party YAML parser with global std specializations,
// so include it outside the namespace wrapper.  The compatibility-patched
// header places only RAD's translator declarations in rrad_rad_core::seqspec.
#include "seqspec.hpp"
