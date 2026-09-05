#pragma once

#include <algorithm>
#include <cstddef>

namespace adapter_thresholds {

constexpr std::size_t kFallbackErrorNumerator = 3;
constexpr std::size_t kFallbackErrorDenominator = 10;

/**
 * Return the shared fallback edit-distance limit for an uncalibrated adapter.
 *
 * Keep this integer-only so layout preparation and read processing cannot
 * drift because of different rounding rules or duplicated constants.
 */
constexpr int fallback_max_edit_distance(std::size_t adapter_length) {
    const auto rounded_up =
        (adapter_length * kFallbackErrorNumerator +
         kFallbackErrorDenominator - 1) /
        kFallbackErrorDenominator;
    return std::max(1, static_cast<int>(rounded_up));
}

static_assert(fallback_max_edit_distance(22) == 7,
              "22 nt adapters must allow seven fallback edits");
static_assert(fallback_max_edit_distance(23) == 7,
              "23 nt adapters must allow seven fallback edits");

}  // namespace adapter_thresholds
