#pragma once

#include <cstdint>

namespace principia {
namespace astronomy {
namespace _solar_system_fingerprints {
namespace internal {

// Indices in the arrays below.
constexpr std::int64_t KSP122 = 0;
constexpr std::int64_t KSP191 = 1;

constexpr std::uint64_t KSPStockSystemFingerprints[] = {
    0x9F3E8BFE0E32C283,
    0x7A077635599E23F2};
constexpr std::uint64_t KSPStabilizedSystemFingerprints[] = {
    0x9F1B6D95399877C6,
    0x4B830CBDF5E77F8D};

}  // namespace internal

using internal::KSP122;
using internal::KSP191;
using internal::KSPStabilizedSystemFingerprints;
using internal::KSPStockSystemFingerprints;

}  // namespace _solar_system_fingerprints
}  // namespace astronomy
}  // namespace principia
