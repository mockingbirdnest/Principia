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

  std::string_view name() const;
  bool IsSet() const;

 private:
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

#define PRINCIPIA_CPUID_FLAG(name, ...)                   \
  inline CPUIDFeatureFlag const name(#name, __VA_ARGS__);

namespace cpuid_feature_flags {
// Table 3-11.
PRINCIPIA_CPUID_FLAG(FPU, 0x01, EDX, 0);    // x87 Floating Point Unit on chip.
PRINCIPIA_CPUID_FLAG(PSN, 0x01, EDX, 18);   // Processor Serial Number.
PRINCIPIA_CPUID_FLAG(SSE, 0x01, EDX, 25);   // Streaming SIMD Extensions.
PRINCIPIA_CPUID_FLAG(SSE2, 0x01, EDX, 26);  // Streaming SIMD Extensions 2.
PRINCIPIA_CPUID_FLAG(SSE3, 0x01, ECX, 0);   // Streaming SIMD Extensions 3.
PRINCIPIA_CPUID_FLAG(FMA, 0x01, ECX, 12);   // Fused Multiply Add.
PRINCIPIA_CPUID_FLAG(SSE4_1, 0x01, ECX, 19);  // Streaming SIMD Extensions 4.1.
PRINCIPIA_CPUID_FLAG(AVX, 0x01, ECX, 28);     // Advanced Vector eXtensions.
// Table 3-8.
PRINCIPIA_CPUID_FLAG(AVX2, 0x07, EBX, 5);       // Advanced Vector eXtensions 2.
PRINCIPIA_CPUID_FLAG(AVX512F, 0x07, EBX, 16);   // AVX-512 Foundation.
PRINCIPIA_CPUID_FLAG(AVX512DQ, 0x07, EBX, 17);  // AVX-512 .
PRINCIPIA_CPUID_FLAG(AVX512_IFMA, 0x07, EBX, 21);  // AVX-512 .
PRINCIPIA_CPUID_FLAG(AVX512PF, 0x07, EBX, 26);     // AVX-512 .
PRINCIPIA_CPUID_FLAG(AVX512ER, 0x07, EBX, 27);     // AVX-512 .
PRINCIPIA_CPUID_FLAG(AVX512CD, 0x07, EBX, 28);     // AVX-512 .
PRINCIPIA_CPUID_FLAG(AVX512BW, 0x07, EBX, 30);     // AVX-512 .
PRINCIPIA_CPUID_FLAG(AVX512VL, 0x07, EBX, 31);     // AVX-512 .
}  // namespace cpuid_feature_flags
}  // namespace internal

using internal::CPUVendorIdentificationString;
namespace cpuid_feature_flags = internal::cpuid_feature_flags

}  // namespace _cpuid
}  // namespace base
}  // namespace principia
