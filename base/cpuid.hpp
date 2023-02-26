#pragma once

#include <cstdint>
#include <string>

namespace principia {
namespace base {
namespace _cpuid {
namespace internal {

// See the Intel® 64 and IA-32 Architectures Software Developer’s Manual,
// Volume 2A, CPUID—CPU Identification.

// Leaf 0.
std::string CPUVendorIdentificationString();

// Leaf 1.
// We represent feature flags as EDX + 2³² ECX.
constexpr std::uint64_t edx_bit = 1;
constexpr std::uint64_t ecx_bit = edx_bit << 32;
enum class CPUFeatureFlags : std::uint64_t {
  // Table 3-11.
  FPU = edx_bit << 0,    // x87 Floating Point Unit on chip.
  PSN = edx_bit << 18,   // Processor Serial Number.
  SSE = edx_bit << 25,   // Streaming SIMD Extensions.
  SSE2 = edx_bit << 26,  // Streaming SIMD Extensions 2.
  // Table 3-10.
  SSE3 = ecx_bit << 0,     // Streaming SIMD Extensions 3.
  FMA = ecx_bit << 12,     // Fused Multiply Add.
  SSE4_1 = ecx_bit << 19,  // Streaming SIMD Extensions 4.1.
  AVX = ecx_bit << 28,     // Advanced Vector eXtensions.
};

// Bitwise or of feature flags; the result represents the union of all features
// in |left| and |right|.
CPUFeatureFlags operator|(CPUFeatureFlags left, CPUFeatureFlags right);

// Whether the CPU has all features listed in |flags|.
bool HasCPUFeatures(CPUFeatureFlags flags);

}  // namespace internal

using internal::CPUFeatureFlags;
using internal::CPUVendorIdentificationString;
using internal::HasCPUFeatures;

}  // namespace _cpuid
}  // namespace base
}  // namespace principia

namespace principia::base {
using namespace principia::base::_cpuid;
}  // namespace principia::base
