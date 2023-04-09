#include "base/cpuid.hpp"

#include <string>

#include "base/macros.hpp"
#if PRINCIPIA_COMPILER_MSVC
#include <intrin.h>
#else
#include <cpuid.h>
#endif

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

std::string CPUVendorIdentificationString() {
  auto const leaf_0 = CPUID(0, 0);
  std::string result(12, '\0');
  std::memcpy(&result[0], &leaf_0.ebx, 4);
  std::memcpy(&result[4], &leaf_0.edx, 4);
  std::memcpy(&result[8], &leaf_0.ecx, 4);
  return result;
}

CPUFeatureFlags operator|(CPUFeatureFlags const left,
                          CPUFeatureFlags const right) {
  return static_cast<CPUFeatureFlags>(static_cast<std::uint64_t>(left) |
                                      static_cast<std::uint64_t>(right));
}

bool HasCPUFeatures(CPUFeatureFlags const flags) {
  auto const leaf_1 = CPUID(1, 0);
  return static_cast<CPUFeatureFlags>(
             (ecx_bit * leaf_1.ecx | edx_bit * leaf_1.edx) &
             static_cast<std::uint64_t>(flags)) == flags;
}

}  // namespace internal
}  // namespace _cpuid
}  // namespace base
}  // namespace principia
