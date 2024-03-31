#include "base/cpuid.hpp"

#include <ranges>
#include <string>
#include <vector>

#include "absl/strings/str_join.h"
#include "base/macros.hpp"  // 🧙 For PRINCIPIA_COMPILER_MSVC.
#include "base/not_null.hpp"
#include "glog/logging.h"

#if PRINCIPIA_COMPILER_MSVC
#include <intrin.h>
#else
#include <cpuid.h>
#endif


namespace principia {
namespace base {
namespace _cpuid {
namespace internal {

using namespace principia::base::_not_null;

namespace {

// This vector is not a static member variable of CPUIDFeatureFlag because we do
// not want to include <vector> in the header.
std::vector<not_null<CPUIDFeatureFlag const*>>& CPUIDFlags() {
  static std::vector<not_null<CPUIDFeatureFlag const*>> result;
  return result;
}

CPUIDResult CPUID(std::uint32_t const eax, std::uint32_t const ecx) {
#if PRINCIPIA_COMPILER_MSVC
  int result[4];
  __cpuidex(result, eax, ecx);
  return {.eax = static_cast<std::uint32_t>(result[0]),
          .ebx = static_cast<std::uint32_t>(result[1]),
          .ecx = static_cast<std::uint32_t>(result[2]),
          .edx = static_cast<std::uint32_t>(result[3])};
#else
  CPUIDResult result;
  __cpuid_count(eax, ecx, result.eax, result.ebx, result.ecx, result.edx);
  return result;
#endif
}

}  // namespace

std::string_view CPUIDFeatureFlag::name() const {
  return name_;
}

bool CPUIDFeatureFlag::IsSet() const {
  return CPUID(leaf_, sub_leaf_).*field_ & (1 << bit_);
}

#define PRINCIPIA_CPUID_FLAG(name, ...) \
  CPUIDFeatureFlag const CPUIDFeatureFlag::name(#name, __VA_ARGS__);

// Table 3-11.
PRINCIPIA_CPUID_FLAG(FPU, 0x01, EDX, 0);  // x87 Floating Point Unit on chip.

PRINCIPIA_CPUID_FLAG(PSN, 0x01, EDX, 18);     // Processor Serial Number.
PRINCIPIA_CPUID_FLAG(SSE, 0x01, EDX, 25);     // Streaming SIMD Extensions.
PRINCIPIA_CPUID_FLAG(SSE2, 0x01, EDX, 26);    // Streaming SIMD Extensions 2.
PRINCIPIA_CPUID_FLAG(SSE3, 0x01, ECX, 0);     // Streaming SIMD Extensions 3.
PRINCIPIA_CPUID_FLAG(FMA, 0x01, ECX, 12);     // Fused Multiply Add.
PRINCIPIA_CPUID_FLAG(SSE4_1, 0x01, ECX, 19);  // Streaming SIMD Extensions 4.1.
PRINCIPIA_CPUID_FLAG(AVX, 0x01, ECX, 28);     // Advanced Vector eXtensions.

// Table 3-8, Structured Extended Feature Flags Enumeration Leaf.
PRINCIPIA_CPUID_FLAG(AVX2, 0x07, EBX, 5);       // Advanced Vector eXtensions 2.
PRINCIPIA_CPUID_FLAG(AVX512F, 0x07, EBX, 16);   // AVX-512 Foundation.
PRINCIPIA_CPUID_FLAG(AVX512DQ, 0x07, EBX, 17);  // DWORD and QWORD instructions.
PRINCIPIA_CPUID_FLAG(AVX512VL, 0x07, EBX, 31);  // Vector Length Extensions.
PRINCIPIA_CPUID_FLAG(AVX512_FP16, 0x07, ECX, 23);  // IEEE-754 binary16.

CPUIDFeatureFlag::CPUIDFeatureFlag(std::string_view const name,
                                   std::uint32_t const leaf,
                                   std::uint32_t const sub_leaf,
                                   std::uint32_t CPUIDResult::*const field,
                                   std::int8_t const bit)
    : name_(name), leaf_(leaf), sub_leaf_(sub_leaf), field_(field), bit_(bit) {
  CHECK_GE(bit, 0);
  CHECK_LT(bit, 32);
  CPUIDFlags().push_back(this);
}

std::string CPUVendorIdentificationString() {
  auto const leaf_0 = CPUID(0, 0);
  std::string result(12, '\0');
  std::memcpy(&result[0], &leaf_0.ebx, 4);
  std::memcpy(&result[4], &leaf_0.edx, 4);
  std::memcpy(&result[8], &leaf_0.ecx, 4);
  return result;
}

std::string ProcessorBrandString() {
  std::string result(48, '\0');
  for (int n = 0; n < 3; ++n) {
    auto const piece = CPUID(0x80000002 + n, 0);
    std::memcpy(&result[n * 16], &piece.eax, 4);
    std::memcpy(&result[n * 16 + 4], &piece.ebx, 4);
    std::memcpy(&result[n * 16 + 8], &piece.ecx, 4);
    std::memcpy(&result[n * 16 + 12], &piece.edx, 4);
  }
  return result;
}

std::string CPUFeatures() {
#if PRINCIPIA_COMPILER_MSVC
  return CPUIDFlags()
      | std::views::filter(&CPUIDFeatureFlag::IsSet)
      | std::views::transform(&CPUIDFeatureFlag::name)  // NOLINT
      | std::views::join_with(' ')
      | std::ranges::to<std::string>();
#else  // TODO(egg): Get rid of this once clang really has C++23.
  std::vector<std::string_view> set_flags;
  for (not_null const flag : CPUIDFlags()) {
    if (flag->IsSet()) {
      set_flags.push_back(flag->name());
    }
  }
  return absl::StrJoin(set_flags, " ");
#endif
}

}  // namespace internal
}  // namespace _cpuid
}  // namespace base
}  // namespace principia
