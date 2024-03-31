#pragma once

#include <cstdint>
#include <string>
#include <string_view>

namespace principia {
namespace base {
namespace _cpuid {
namespace internal {

struct CPUIDResult {
  std::uint32_t eax;
  std::uint32_t ebx;
  std::uint32_t ecx;
  std::uint32_t edx;
};

class CPUIDFeatureFlag {
 public:
  CPUIDFeatureFlag(CPUIDFeatureFlag const&) = delete;

  std::string_view name() const;
  bool IsSet() const;

  static const CPUIDFeatureFlag FPU;       // x87 Floating Point Unit on chip.
  static const CPUIDFeatureFlag PSN;       // Processor Serial Number.
  static const CPUIDFeatureFlag SSE;       // Streaming SIMD Extensions.
  static const CPUIDFeatureFlag SSE2;      // Streaming SIMD Extensions 2.
  static const CPUIDFeatureFlag SSE3;      // Streaming SIMD Extensions 3.

  static const CPUIDFeatureFlag FMA;          // Fused Multiply Add.
  static const CPUIDFeatureFlag SSE4_1;       // Streaming SIMD Extensions 4.1.
  static const CPUIDFeatureFlag AVX;          // Advanced Vector eXtensions.
  static const CPUIDFeatureFlag AVX2;         // Advanced Vector eXtensions 2.
  static const CPUIDFeatureFlag AVX512F;      // AVX-512 Foundation.
  static const CPUIDFeatureFlag AVX512DQ;     // DWORD and QWORD instructions.
  static const CPUIDFeatureFlag AVX512VL;     // Vector Length Extensions.
  static const CPUIDFeatureFlag AVX512_FP16;  // IEEE-754 binary16.

 private:
  CPUIDFeatureFlag(std::string_view const name,
                   std::uint32_t const leaf,
                   std::uint32_t const sub_leaf,
                   std::uint32_t CPUIDResult::*const field,
                   std::int8_t const bit);
  CPUIDFeatureFlag(std::string_view const name,
                   std::uint32_t const leaf,
                   std::uint32_t CPUIDResult::*const field,
                   std::int8_t const bit)
      : CPUIDFeatureFlag(name, leaf, 0, field, bit) {}

  std::string_view name_;
  std::uint32_t leaf_;
  std::uint32_t sub_leaf_;
  std::uint32_t CPUIDResult::*field_;
  std::int8_t bit_;
};

constexpr std::uint32_t CPUIDResult::*EAX = &CPUIDResult::eax;
constexpr std::uint32_t CPUIDResult::*EBX = &CPUIDResult::ebx;
constexpr std::uint32_t CPUIDResult::*ECX = &CPUIDResult::ecx;
constexpr std::uint32_t CPUIDResult::*EDX = &CPUIDResult::edx;

std::string CPUVendorIdentificationString();
std::string ProcessorBrandString();

std::string CPUFeatures();

}  // namespace internal

#undef PRINCIPIA_CPUID_FLAG

using internal::CPUFeatures;
using internal::CPUIDFeatureFlag;
using internal::CPUVendorIdentificationString;
using internal::ProcessorBrandString;

}  // namespace _cpuid
}  // namespace base
}  // namespace principia
