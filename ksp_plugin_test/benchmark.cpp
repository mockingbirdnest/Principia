
#include "ksp_plugin/plugin.hpp"

#include <string>
#include <vector>

#include "base/push_deserializer.hpp"
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
using base::PushDeserializer;
using geometry::Instant;
using interface::principia__AdvanceTime;
using interface::principia__DeletePlugin;
using interface::principia__DeserializePluginHexadecimal;
using interface::principia__FutureCatchUpVessel;
using interface::principia__FutureWaitForVesselToCatchUp;
using interface::principia__IteratorDelete;
using quantities::Frequency;
using quantities::Time;
using quantities::si::Hertz;
using quantities::si::Second;
using testing_utilities::ReadFromBinaryFile;
using testing_utilities::ReadFromHexadecimalFile;

namespace ksp_plugin {

void BM_PluginIntegrationBenchmark(benchmark::State& state) {
  auto const plugin = Plugin::ReadFromMessage(
      ParseFromBytes<serialization::Plugin>(ReadFromBinaryFile(
          SOLUTION_DIR / "ksp_plugin_test" / "3 vessels.proto.bin")));

  std::vector<GUID> const vessel_guids = {
      "70ff8dc0-a4dd-4b8c-868b-35ddb01e32bc",
      "abd95a7e-6b8b-4dba-a1e9-c96cd594cd67",
      "b86d2efd-5150-4a44-8c36-04820a85e861"};

  static constexpr int warp_factor = 6E6;
  static constexpr Frequency refresh_frequency = 50 * Hertz;
  static constexpr Time step = warp_factor / refresh_frequency;
  for (auto _ : state) {
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

void BM_PluginDeserializationBenchmark(benchmark::State& state) {
  char const compressor[] = "";
  std::string const hexadecimal_simple_plugin(
      ReadFromHexadecimalFile(
          SOLUTION_DIR / "ksp_plugin_test" / "simple_plugin.proto.hex"));

  int bytes_processed = 0;
  for (auto _ : state) {
    PushDeserializer* deserializer = nullptr;
    Plugin const* plugin = nullptr;
    principia__DeserializePluginHexadecimal(hexadecimal_simple_plugin.c_str(),
                                            hexadecimal_simple_plugin.size(),
                                            &deserializer,
                                            &plugin,
                                            compressor);
    principia__DeserializePluginHexadecimal(hexadecimal_simple_plugin.c_str(),
                                            0,
                                            &deserializer,
                                            &plugin,
                                            compressor);
    principia__DeletePlugin(&plugin);
    bytes_processed += hexadecimal_simple_plugin.size() >> 1;
  }
  state.SetBytesProcessed(bytes_processed);
}

BENCHMARK(BM_PluginDeserializationBenchmark);
BENCHMARK(BM_PluginIntegrationBenchmark);

// .\Release\x64\ksp_plugin_test_tests.exe --gtest_filter=PluginBenchmark.DISABLED_All --gtest_also_run_disabled_tests  // NOLINT
TEST(PluginBenchmark, DISABLED_All) {
  benchmark::RunSpecifiedBenchmarks();
}

}  // namespace ksp_plugin
}  // namespace principia
