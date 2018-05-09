
#include "ksp_plugin/plugin.hpp"

#include <string>
#include <vector>

#include "base/serialization.hpp"
#include "benchmark/benchmark.h"
#include "geometry/named_quantities.hpp"
#include "gtest/gtest.h"
#include "ksp_plugin/interface.hpp"
#include "quantities/quantities.hpp"
#include "serialization/ksp_plugin.pb.h"
#include "testing_utilities/serialization.hpp"

namespace principia {

using base::ParseFromBytes;
using geometry::Instant;
using interface::principia__AdvanceTime;
using interface::principia__FutureCatchUpVessel;
using interface::principia__FutureWaitForVesselToCatchUp;
using interface::principia__IteratorDelete;
using quantities::Frequency;
using quantities::Time;
using quantities::si::Hertz;
using quantities::si::Second;
using testing_utilities::ReadFromBinaryFile;

namespace ksp_plugin {

void BM_PluginIntegrationBenchmark(benchmark::State& state) {
  auto const binary_plugin = ReadFromBinaryFile(
      SOLUTION_DIR / "ksp_plugin_test" / "3 vessels.proto.bin");
  serialization::Plugin serialized_plugin;
  // TODO(phl): For some reason this doesn't CHECK because it wants to read
  // beyond the end |binary_plugin|.  Figure out why (could it be an extra
  // byte?).
  serialized_plugin.ParseFromArray(binary_plugin.data(), binary_plugin.size());
  auto const plugin = Plugin::ReadFromMessage(serialized_plugin);

  std::vector<GUID> const vessel_guids = {
      "70ff8dc0-a4dd-4b8c-868b-35ddb01e32bc",
      "abd95a7e-6b8b-4dba-a1e9-c96cd594cd67",
      "b86d2efd-5150-4a44-8c36-04820a85e861"};

  static constexpr int warp_factor = 6E6;
  static constexpr Frequency refresh_frequency = 50 * Hertz;
  static constexpr Time step = warp_factor / refresh_frequency;
  while (state.KeepRunning()) {
    principia__AdvanceTime(
        plugin.get(),
        (plugin->CurrentTime() + step - plugin->GameEpoch()) / Second,
        /*planetarium_rotation=*/45);
    std::vector<PileUpFuture*> futures;
    for (GUID const& vessel_guid : vessel_guids) {
      futures.push_back(
          principia__FutureCatchUpVessel(plugin.get(), vessel_guid.c_str()));
    }
    for (auto& future : futures) {
      Iterator* iterator;
      principia__FutureWaitForVesselToCatchUp(plugin.get(),
                                              &future,
                                              &iterator);
      principia__IteratorDelete(&iterator);
    }
  }
}

BENCHMARK(BM_PluginIntegrationBenchmark);

TEST(PluginBenchmark, DISABLED_3Vessels) {
  benchmark::RunSpecifiedBenchmarks();
}

}  // namespace ksp_plugin
}  // namespace principia
