#include "nanobenchmarks/microarchitectures.hpp"

#include <map>
#include <regex>
#include <utility>
#include <vector>

#include "base/cpuid.hpp"
#include "base/macros.hpp"  // 🧙 For PRINCIPIA_COMPILER_CLANG.
#include "glog/logging.h"

BENCHMARK_EXTERN_C_FUNCTION(identity);
BENCHMARK_EXTERN_C_FUNCTION(sqrtps_xmm0_xmm0);
BENCHMARK_EXTERN_C_FUNCTION(sqrtsd_xmm0_xmm0);
BENCHMARK_EXTERN_C_FUNCTION(mulsd_xmm0_xmm0);
BENCHMARK_EXTERN_C_FUNCTION(mulsd_xmm0_xmm0_4x);
#if PRINCIPIA_COMPILER_CLANG
asm(R"(
.intel_syntax
identity:
  ret
sqrtps_xmm0_xmm0:
  sqrtps xmm0, xmm0
  ret
sqrtsd_xmm0_xmm0:
  sqrtsd xmm0, xmm0
  ret
mulsd_xmm0_xmm0:
  mulsd xmm0, xmm0
  ret
mulsd_xmm0_xmm0_4x:
  mulsd xmm0, xmm0
  mulsd xmm0, xmm0
  mulsd xmm0, xmm0
  mulsd xmm0, xmm0
  ret
)");
#endif

namespace principia {
namespace nanobenchmarks {
namespace _microarchitectures {
namespace internal {

using namespace principia::base::_cpuid;

namespace {
static std::vector<
    std::pair<std::regex, std::map<BenchmarkedFunction, int>>> const&
    microarchitectures = *new std::vector{
        // Skylake, Cascade Lake, Coffee Lake, Cannon Lake, Ice Lake, Tiger
        // Lake, Golden Cove(?).
        std::pair{std::regex(R"(((6|7|9|10|11|12)th Gen Intel\(R\) Core\(TM\))"
                             R"(|Intel\(R\) Xeon\(R\) W-[23]).*)"),
                  std::map{std::pair{&identity, 0},
                           std::pair{&mulsd_xmm0_xmm0, 4},
                           std::pair{&mulsd_xmm0_xmm0_4x, 4 * 4},
                           std::pair{&sqrtps_xmm0_xmm0, 12}}},
        // Zen3.
        std::pair{std::regex("AMD Ryzen Threadripper PRO 5.*"),
                  std::map{std::pair{&identity, 0},
                           std::pair{&mulsd_xmm0_xmm0, 3},
                           std::pair{&mulsd_xmm0_xmm0_4x, 4 * 3},
                           std::pair{&sqrtps_xmm0_xmm0, 14}}}};
}  // namespace

std::map<BenchmarkedFunction, int> const& ReferenceCycleCounts() {
  static std::map<BenchmarkedFunction, int> const& result = [] {
    for (auto const& [regex, architecture] : microarchitectures) {
      if (std::regex_match(ProcessorBrandString(), regex)) {
        return architecture;
      }
    }
    LOG(FATAL) << "Unknown architecture " << CPUVendorIdentificationString()
               << " " << ProcessorBrandString();
  }();
  return result;
}

}  // namespace internal
}  // namespace _microarchitectures
}  // namespace nanobenchmarks
}  // namespace principia
