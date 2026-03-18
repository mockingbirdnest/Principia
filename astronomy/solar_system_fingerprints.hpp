#pragma once

#include <array>
#include <cstdint>

namespace principia {
namespace astronomy {
namespace _solar_system_fingerprints {
namespace internal {

// Indices in the arrays below.
constexpr std::int64_t KSP122 = 0;
constexpr std::int64_t KSP191PreLegendre = 1;
constexpr std::int64_t KSP191 = 2;

constexpr std::array KSPStockSystemFingerprints{
    0x9F3E8BFE0E32C283ull,
    0x7A077635599E23F2ull,
    0xFE67F3BAEE725803ull};
constexpr std::array KSPStabilizedSystemFingerprints{
    0x9F1B6D95399877C6ull,
    0x4B830CBDF5E77F8Dull,
    0xB1BA690A45CAD577ull};

}  // namespace internal

using internal::KSP122;
using internal::KSP191;
using internal::KSP191PreLegendre;
using internal::KSPStabilizedSystemFingerprints;
using internal::KSPStockSystemFingerprints;

}  // namespace _solar_system_fingerprints
}  // namespace astronomy
}  // namespace principia
