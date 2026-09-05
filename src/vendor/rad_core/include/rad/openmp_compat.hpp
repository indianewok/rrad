#pragma once

// OpenMP is optional for the embedded backend. In particular, CRAN's macOS
// toolchain may not provide an omp.h/runtime pair compatible with the system
// Apple clang. Keep the call sites usable in a deterministic serial build
// without requiring an external compiler runtime.
#ifdef _OPENMP
#include <omp.h>
#else
inline int omp_get_max_threads() noexcept { return 1; }
inline int omp_get_num_threads() noexcept { return 1; }
inline int omp_get_thread_num() noexcept { return 0; }
inline void omp_set_num_threads(int) noexcept {}
#endif
