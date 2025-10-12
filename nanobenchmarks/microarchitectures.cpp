#include "nanobenchmarks/microarchitectures.hpp"

#include <memory>
#include <regex>
#include <string>
#include <utility>
#include <vector>

#include "absl/container/btree_map.h"
#include "absl/container/flat_hash_map.h"
#include "base/cpuid.hpp"
#include "base/macros.hpp"  // 🧙 For PRINCIPIA_COMPILER_CLANG.
#include "glog/logging.h"
#include "nanobenchmarks/nanobenchmark.hpp"

namespace principia {
namespace nanobenchmarks {
namespace _microarchitectures {
namespace internal {

using namespace principia::base::_cpuid;
using namespace principia::nanobenchmarks::_nanobenchmark;

NANOBENCHMARK_EXTERN_C_FUNCTION(maxps_xmm0_xmm0);
NANOBENCHMARK_EXTERN_C_FUNCTION(sqrtps_xmm0_xmm0);
NANOBENCHMARK_EXTERN_C_FUNCTION(sqrtsd_xmm0_xmm0);
NANOBENCHMARK_EXTERN_C_FUNCTION(mulsd_xmm0_xmm0);
NANOBENCHMARK_EXTERN_C_FUNCTION(mulsd_xmm0_xmm0_4x);

#if PRINCIPIA_COMPILER_CLANG
asm(R"(
.intel_syntax
_maxps_xmm0_xmm0:
  maxps xmm0, xmm0
  ret
_sqrtps_xmm0_xmm0:
  sqrtps xmm0, xmm0
  ret
_sqrtsd_xmm0_xmm0:
  sqrtsd xmm0, xmm0
  ret
_mulsd_xmm0_xmm0:
  mulsd xmm0, xmm0
  ret
_mulsd_xmm0_xmm0_4x:
  mulsd xmm0, xmm0
  mulsd xmm0, xmm0
  mulsd xmm0, xmm0
  mulsd xmm0, xmm0
  ret
)");
#endif

namespace {
static std::vector<
    std::pair<std::regex, absl::flat_hash_map<BenchmarkedFunction, int>>> const&
    microarchitectures = *new std::vector{
        // Skylake, Cascade Lake, Coffee Lake, Cannon Lake, Ice Lake, Tiger
        // Lake, Golden Cove(?).
        std::pair{std::regex(R"(((6|7|9|10|11|12)th Gen Intel\(R\) Core\(TM\))"
                             R"(|Intel\(R\) Xeon\(R\) W-[23]).*)"),
                  absl::flat_hash_map<BenchmarkedFunction, int>{
                      std::pair{&maxps_xmm0_xmm0, 1},
                      std::pair{&mulsd_xmm0_xmm0, 4},
                      std::pair{&mulsd_xmm0_xmm0_4x, 4 * 4},
                      std::pair{&sqrtps_xmm0_xmm0, 12}}},
        // Zen3.
        std::pair{std::regex("AMD Ryzen Threadripper PRO 5.*"),
                  absl::flat_hash_map<BenchmarkedFunction, int>{
                      std::pair{&maxps_xmm0_xmm0, 1},
                      std::pair{&mulsd_xmm0_xmm0, 3},
                      std::pair{&mulsd_xmm0_xmm0_4x, 4 * 3},
                      std::pair{&sqrtps_xmm0_xmm0, 14}}},
        // Rosetta 2.
        std::pair{std::regex("VirtualApple .*"),
                  absl::flat_hash_map<BenchmarkedFunction, int>{
                      std::pair{&maxps_xmm0_xmm0, 1},
                      std::pair{&mulsd_xmm0_xmm0, 4},
                      std::pair{&mulsd_xmm0_xmm0_4x, 4 * 4}}}};
}  // namespace

std::vector<NanobenchmarkAndCycles> const& ReferenceCycleCounts() {
  static std::vector<NanobenchmarkAndCycles>* const reference_cycle_counts =
      [] {
        for (auto const& [regex, architecture] : microarchitectures) {
          if (std::regex_match(ProcessorBrandString(), regex)) {
            absl::btree_map<std::string, NanobenchmarkAndCycles>
                sorted_reference_cycle_counts;
            for (auto const& [function, cycles] : architecture) {
              auto const* const nanobenchmark =
                  NanobenchmarkRegistry::NanobenchmarkFor(function);
              CHECK(nanobenchmark != nullptr)
                  << "No nanobenchmark for function at " << function;
              sorted_reference_cycle_counts.emplace(
                  nanobenchmark->name(),
                  NanobenchmarkAndCycles{nanobenchmark, cycles});
            }
            return new std::vector<NanobenchmarkAndCycles>(
                std::from_range,
                sorted_reference_cycle_counts | std::views::values);
          }
        }
        LOG(FATAL) << "Unknown architecture " << CPUVendorIdentificationString()
                   << " " << ProcessorBrandString();
      }();
  return *reference_cycle_counts;
}

}  // namespace internal
}  // namespace _microarchitectures
}  // namespace nanobenchmarks
}  // namespace principia
