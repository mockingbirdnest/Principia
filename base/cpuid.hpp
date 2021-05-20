#pragma once

#include <cstdint>

#include <string>

namespace principia {
namespace base {
namespace internal_cpuid {

// See the Intel® 64 and IA-32 Architectures Software Developer’s Manual,
// Volume 2A, CPUID—CPU Identification.

// Leaf 0.
std::string VendorIdentificationString();

// Leaf 1.
// We represent feature flags as EDX + 2³² ECX.
constexpr std::uint64_t edx_bit = 1;
constexpr std::uint64_t ecx_bit = edx_bit << 32;
enum class FeatureFlags : std::uint64_t {
  // Table 3-11.
  FPU = edx_bit << 0,
  SSE = edx_bit << 25,
  SSE2 = edx_bit << 26,
  // Table 3-10.
  SSE3 = ecx_bit << 0,
  FMA = ecx_bit << 12,
  SSE4_1 = ecx_bit << 19,
  AVX = ecx_bit << 28,
  NotUsed = ecx_bit << 31,  // Always 0.
};

FeatureFlags operator|(FeatureFlags left, FeatureFlags right);

bool HasCPUFeatures(FeatureFlags const flags);

}  // namespace internal_cpuid

using internal_cpuid::VendorIdentificationString;
using internal_cpuid::FeatureFlags;
using internal_cpuid::HasCPUFeatures;

}  // namespace base
}  // namespace principia
