#pragma once

#include "nanobenchmarks/microarchitectures.hpp"

#include <memory>
#include <regex>
#include <string>
#include <utility>
#include <vector>

#include "absl/container/btree_map.h"
#include "base/cpuid.hpp"
#include "base/macros.hpp"  // 🧙 For PRINCIPIA_COMPILER_CLANG.
#include "geometry/instant.hpp"
#include "geometry/space.hpp"
#include "glog/logging.h"
#include "nanobenchmarks/dependencies.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace nanobenchmarks {
namespace _microarchitectures {
namespace internal {

using namespace principia::base::_cpuid;
using namespace principia::geometry::_instant;
using namespace principia::geometry::_space;
using namespace principia::nanobenchmarks::_dependencies;
using namespace principia::physics::_degrees_of_freedom;
using namespace principia::quantities::_si;

NANOBENCHMARK_EXTERN_C_FUNCTION(identity);
NANOBENCHMARK_EXTERN_C_FUNCTION(sqrtps_xmm0_xmm0);
NANOBENCHMARK_EXTERN_C_FUNCTION(sqrtsd_xmm0_xmm0);
NANOBENCHMARK_EXTERN_C_FUNCTION(mulsd_xmm0_xmm0);
NANOBENCHMARK_EXTERN_C_FUNCTION(mulsd_xmm0_xmm0_4x);

using DisplacementInstantNanobenchmark =
    Nanobenchmark<Displacement<World>, Instant>;
using RelativeDegreesOfFreedomInstantNanobenchmark =
    Nanobenchmark<RelativeDegreesOfFreedom<World>, Instant>;

#if 0
// NOTE: To figure out the calling convention or to get the assembly code for a
// simple nanobenchmark, uncomment this code and do the necessary adjustments,
// and then look in `main.asm` for the `PROC` for
// `NanobenchmarkCase@RelativeDegreesOfFreedomInstantNanobenchmark_Reference_Nanobenchmark`.
NANOBENCHMARK_FIXTURE(RelativeDegreesOfFreedomInstantNanobenchmark, Reference) {
  double const x = (argument - Instant()) / Second;
  return RelativeDegreesOfFreedom<World>(
      Displacement<World>({x * Metre, x * Metre, x * Metre}),
      Velocity<World>(
          {x * Metre / Second, x * Metre / Second, x * Metre / Second}));
}
#endif

// NOTE: To find the mangled name, just add a macro with a bogus full name, and
// copy the mangled name from the error message.
NANOBENCHMARK_EXTERN_ALTERNATE_NAME_FUNCTION(
    DisplacementInstantNanobenchmark,
    fill3_from_xmm0,
    "?fill3_from_xmm0@internal@_microarchitectures@nanobenchmarks@principia@@"
    "YQ?AV?$Multivector@V?$Quantity@U?$Dimensions@$00$0A@$0A@$0A@$0A@$0A@$0A@$"
    "0A@@internal@_dimensions@quantities@principia@@@internal@_quantities@"
    "quantities@principia@@U?$Frame@W4Frame_TestTag@serialization@principia@@$"
    "0A@$00$00@2_frame@geometry@5@$00@1_grassmann@geometry@4@V?$Point@V?$"
    "Quantity@U?$Dimensions@$0A@$0A@$00$0A@$0A@$0A@$0A@$0A@@internal@_"
    "dimensions@quantities@principia@@@internal@_quantities@quantities@"
    "principia@@@1_point@74@@Z");
NANOBENCHMARK_EXTERN_ALTERNATE_NAME_FUNCTION(
    DisplacementInstantNanobenchmark,
    fill_from_ymm0,
    "?fill_from_ymm0@internal@_microarchitectures@nanobenchmarks@principia@@YQ?"
    "AV?$Multivector@V?$Quantity@U?$Dimensions@$00$0A@$0A@$0A@$0A@$0A@$0A@$0A@@"
    "internal@_dimensions@quantities@principia@@@internal@_quantities@"
    "quantities@principia@@U?$Frame@W4Frame_TestTag@serialization@principia@@$"
    "0A@$00$00@2_frame@geometry@5@$00@1_grassmann@geometry@4@V?$Point@V?$"
    "Quantity@U?$Dimensions@$0A@$0A@$00$0A@$0A@$0A@$0A@$0A@@internal@_"
    "dimensions@quantities@principia@@@internal@_quantities@quantities@"
    "principia@@@1_point@74@@Z");
NANOBENCHMARK_EXTERN_ALTERNATE_NAME_FUNCTION(
    RelativeDegreesOfFreedomInstantNanobenchmark,
    fill6_from_xmm0,
    "?fill6_from_xmm0@internal@_microarchitectures@nanobenchmarks@principia@@"
    "YQ?AV?$Pair@V?$Multivector@V?$Quantity@U?$Dimensions@$00$0A@$0A@$0A@$0A@$"
    "0A@$0A@$0A@@internal@_dimensions@quantities@principia@@@internal@_"
    "quantities@quantities@principia@@U?$Frame@W4Frame_TestTag@serialization@"
    "principia@@$0A@$00$00@2_frame@geometry@5@$00@internal@_grassmann@geometry@"
    "principia@@V?$Multivector@V?$Quantity@U?$Dimensions@$00$0A@$0?0$0A@$0A@$"
    "0A@$0A@$0A@@internal@_dimensions@quantities@principia@@@internal@_"
    "quantities@quantities@principia@@U?$Frame@W4Frame_TestTag@serialization@"
    "principia@@$0A@$00$00@2_frame@geometry@5@$00@2345@@1_pair@geometry@4@V?$"
    "Point@V?$Quantity@U?$Dimensions@$0A@$0A@$00$0A@$0A@$0A@$0A@$0A@@internal@_"
    "dimensions@quantities@principia@@@internal@_quantities@quantities@"
    "principia@@@1_point@74@@Z");
NANOBENCHMARK_EXTERN_ALTERNATE_NAME_FUNCTION(
    RelativeDegreesOfFreedomInstantNanobenchmark,
    fill_from_ymm0_ymm1,
    "?fill_from_ymm0_ymm1@internal@_microarchitectures@nanobenchmarks@"
    "principia@@YQ?AV?$Pair@V?$Multivector@V?$Quantity@U?$Dimensions@$00$0A@$"
    "0A@$0A@$0A@$0A@$0A@$0A@@internal@_dimensions@quantities@principia@@@"
    "internal@_quantities@quantities@principia@@U?$Frame@W4Frame_TestTag@"
    "serialization@principia@@$0A@$00$00@2_frame@geometry@5@$00@internal@_"
    "grassmann@geometry@principia@@V?$Multivector@V?$Quantity@U?$Dimensions@$"
    "00$0A@$0?0$0A@$0A@$0A@$0A@$0A@@internal@_dimensions@quantities@principia@@"
    "@internal@_quantities@quantities@principia@@U?$Frame@W4Frame_TestTag@"
    "serialization@principia@@$0A@$00$00@2_frame@geometry@5@$00@2345@@1_pair@"
    "geometry@4@V?$Point@V?$Quantity@U?$Dimensions@$0A@$0A@$00$0A@$0A@$0A@$0A@$"
    "0A@@internal@_dimensions@quantities@principia@@@internal@_quantities@"
    "quantities@principia@@@1_point@74@@Z");

#if PRINCIPIA_COMPILER_CLANG
#if OS_MACOSX
asm(R"(
.intel_syntax
_identity:
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
#elif OS_LINUX
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
_ZN9principia14nanobenchmarks19_microarchitectures8internal15fill3_from_xmm0ENS_8geometry6_point8internal5PointINS_10quantities11_quantities8internal8QuantityINS7_11_dimensions8internal10DimensionsILl0ELl0ELl1ELl0ELl0ELl0ELl0ELl0EEEEEEE:
  vxorpd      xmm1, xmm1, xmm1
  vunpcklpd   xmm2, xmm0, xmm1
  vmovddup    xmm1, xmm2
  vinsertf128 ymm2, ymm1, xmm2, 1
  vmovupd     YMMWORD PTR [rcx], ymm2
  mov         rax, rcx
  vzeroupper
  ret
_ZN9principia14nanobenchmarks19_microarchitectures8internal14fill_from_ymm0ENS_8geometry6_point8internal5PointINS_10quantities11_quantities8internal8QuantityINS7_11_dimensions8internal10DimensionsILl0ELl0ELl1ELl0ELl0ELl0ELl0ELl0EEEEEEE:
  vmovupd     YMMWORD PTR [rcx], ymm0
  mov         rax, rcx
  vzeroupper
  ret
_ZN9principia14nanobenchmarks19_microarchitectures8internal15fill6_from_xmm0ENS_8geometry6_point8internal5PointINS_10quantities11_quantities8internal8QuantityINS7_11_dimensions8internal10DimensionsILl0ELl0ELl1ELl0ELl0ELl0ELl0ELl0EEEEEEE:
  vxorpd      xmm1, xmm1, xmm1
  vunpcklpd   xmm2, xmm0, xmm1
  vmovddup    xmm1, xmm2
  vinsertf128 ymm2, ymm1, xmm2, 1
  vmovupd     YMMWORD PTR [rcx], ymm2
  vmovupd     YMMWORD PTR [rcx+32], ymm2
  mov         rax, rcx
  vzeroupper
  ret
_ZN9principia14nanobenchmarks19_microarchitectures8internal19fill_from_ymm0_ymm1ENS_8geometry6_point8internal5PointINS_10quantities11_quantities8internal8QuantityINS7_11_dimensions8internal10DimensionsILl0ELl0ELl1ELl0ELl0ELl0ELl0ELl0EEEEEEE:
  vmovupd     YMMWORD PTR [rcx], ymm0
  vmovupd     YMMWORD PTR [rcx+32], ymm1
  mov         rax, rcx
  vzeroupper
  ret
)");
#endif
#endif

namespace {

// Skylake, Cascade Lake, Coffee Lake, Cannon Lake, Ice Lake, Tiger Lake, Golden
// Cove(?).
std::regex const intel(R"(((6|7|9|10|11|12)th Gen Intel\(R\) Core\(TM\))"
                       R"(|Intel\(R\) Xeon\(R\) W-[23]).*)");
std::regex const zen3("AMD Ryzen Threadripper PRO 5.*");
std::regex const rosetta2("VirtualApple .*");

// Don't use `absl::flat_hash_map` below, for some reason it introduces noise in
// the timings.
template<typename Value, typename Argument>
using MicroarchitectureDescription = std::vector<
    std::pair<std::regex,
              absl::btree_map<
                  typename Nanobenchmark<Value, Argument>::BenchmarkedFunction,
                  int>>>;

template<typename Value, typename Argument>
MicroarchitectureDescription<Value, Argument> const microarchitectures;

template<>
MicroarchitectureDescription<double, double> const
    microarchitectures<double, double>{
        std::pair{
            intel,
            absl::btree_map<Nanobenchmark<double, double>::BenchmarkedFunction,
                            int>{std::pair{&identity, 0},
                                 std::pair{&mulsd_xmm0_xmm0, 4},
                                 std::pair{&mulsd_xmm0_xmm0_4x, 4 * 4},
                                 std::pair{&sqrtps_xmm0_xmm0, 12}}},
        std::pair{
            zen3,
            absl::btree_map<Nanobenchmark<double, double>::BenchmarkedFunction,
                            int>{std::pair{&identity, 0},
                                 std::pair{&mulsd_xmm0_xmm0, 3},
                                 std::pair{&mulsd_xmm0_xmm0_4x, 4 * 3},
                                 std::pair{&sqrtps_xmm0_xmm0, 14}}},
        std::pair{
            rosetta2,
            absl::btree_map<Nanobenchmark<double, double>::BenchmarkedFunction,
                            int>{std::pair{&identity, 0},
                                 std::pair{&mulsd_xmm0_xmm0, 4},
                                 std::pair{&mulsd_xmm0_xmm0_4x, 4 * 4}}}};

template<>
MicroarchitectureDescription<Displacement<World>, Instant> const
    microarchitectures<Displacement<World>, Instant>{
        std::pair{
            intel,
            absl::btree_map<Nanobenchmark<Displacement<World>,
                                          Instant>::BenchmarkedFunction,
                            int>{std::pair{&fill3_from_xmm0, 1 + 1 + 3 + 3},
                                 std::pair{&fill_from_ymm0, 3}}},
        std::pair{
            zen3,
            absl::btree_map<Nanobenchmark<Displacement<World>,
                                          Instant>::BenchmarkedFunction,
                            int>{std::pair{&fill3_from_xmm0, 1 + 1 + 2 + 5},
                                 std::pair{&fill_from_ymm0, 5}}},
        std::pair{rosetta2,
                  absl::btree_map<Nanobenchmark<Displacement<World>,
                                                Instant>::BenchmarkedFunction,
                                  int>{std::pair{&fill3_from_xmm0, 0},
                                       std::pair{&fill_from_ymm0, 0}}}};

template<>
MicroarchitectureDescription<RelativeDegreesOfFreedom<World>, Instant> const
    microarchitectures<RelativeDegreesOfFreedom<World>, Instant>{
        std::pair{
            intel,
            absl::btree_map<Nanobenchmark<RelativeDegreesOfFreedom<World>,
                                          Instant>::BenchmarkedFunction,
                            int>{std::pair{&fill6_from_xmm0, 1 + 1 + 3 + 3},
                                 std::pair{&fill_from_ymm0_ymm1, 3}}},
        std::pair{
            zen3,
            absl::btree_map<Nanobenchmark<RelativeDegreesOfFreedom<World>,
                                          Instant>::BenchmarkedFunction,
                            int>{std::pair{&fill6_from_xmm0, 1 + 1 + 2 + 5},
                                 std::pair{&fill_from_ymm0_ymm1, 5}}},
        std::pair{rosetta2,
                  absl::btree_map<Nanobenchmark<RelativeDegreesOfFreedom<World>,
                                                Instant>::BenchmarkedFunction,
                                  int>{std::pair{&fill6_from_xmm0, 0},
                                       std::pair{&fill_from_ymm0_ymm1, 0}}}};

}  // namespace

template<typename Value, typename Argument>
std::vector<NanobenchmarkAndCycles<Value, Argument>> const&
ReferenceCycleCounts() {
  static std::vector<NanobenchmarkAndCycles<Value, Argument>>* const
      reference_cycle_counts = [] {
        for (auto const& [regex, architecture] :
             microarchitectures<Value, Argument>) {
          if (std::regex_match(ProcessorBrandString(), regex)) {
            absl::btree_map<std::string,
                            NanobenchmarkAndCycles<Value, Argument>>
                sorted_reference_cycle_counts;
            for (auto const& [function, cycles] : architecture) {
              auto const* const nanobenchmark =
                  NanobenchmarkRegistry<Value, Argument>::NanobenchmarkFor(
                      function);
              CHECK(nanobenchmark != nullptr)
                  << "No nanobenchmark for function at " << function;
              sorted_reference_cycle_counts.emplace(
                  nanobenchmark->name(),
                  NanobenchmarkAndCycles{nanobenchmark, cycles});
            }
            return new std::vector<NanobenchmarkAndCycles<Value, Argument>>(
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
