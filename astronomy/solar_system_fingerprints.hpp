#pragma once

#include <cstdint>

namespace principia {
namespace astronomy {

// Indices in the arrays below.
constexpr std::int64_t KSP122 = 0;
constexpr std::int64_t KSP191 = 1;

constexpr std::uint64_t KSPStockSystemFingerprints[] = {
    0xC366DDFC05246F42,
    0x7A077635599E23F2};
constexpr std::uint64_t KSPStabilizedSystemFingerprints[] = {
    0xC6E30693245BE096,
    0x4B830CBDF5E77F8D};

}  // namespace astronomy
}  // namespace principia
