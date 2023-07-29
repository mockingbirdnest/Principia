#include "ksp_plugin/flight_plan_optimizer.hpp"

#include "astronomy/date_time.hpp"
#include "astronomy/time_scales.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "ksp_plugin/plugin.hpp"
#include "ksp_plugin_test/plugin_io.hpp"

namespace principia {
namespace ksp_plugin {

using ::testing::ElementsAre;
using ::testing::Eq;
using ::testing::Pair;
using namespace principia::astronomy::_date_time;
using namespace principia::astronomy::_time_scales;
using namespace principia::base::_not_null;
using namespace principia::ksp_plugin::_flight_plan_optimizer;
using namespace principia::ksp_plugin::_plugin;
using namespace principia::ksp_plugin_test::_plugin_io;
using namespace principia::quantities::_named_quantities;
using namespace principia::quantities::_si;

class FlightPlanOptimizerTest : public testing::Test {};

TEST_F(FlightPlanOptimizerTest, Reach) {
  not_null<std::unique_ptr<Plugin const>> plugin = ReadPluginFromFile(
      SOLUTION_DIR / "ksp_plugin_test" / "saves" / "3072.proto.b64",
      /*compressor=*/"gipfeli",
      /*decoder=*/"base64");

  auto const ifnity = plugin->GetVessel("29142a79-7acd-47a9-a34d-f9f2a8e1b4ed");
  EXPECT_THAT(ifnity->name(), Eq("IFNITY-5.2"));
  EXPECT_THAT(TTSecond(ifnity->trajectory().front().time),
              Eq("1970-08-14T08:03:46"_DateTime));
  EXPECT_THAT(TTSecond(ifnity->psychohistory()->back().time),
              Eq("1970-08-14T08:47:05"_DateTime));
  ASSERT_TRUE(ifnity->has_flight_plan());
  ifnity->ReadFlightPlanFromMessage();
  EXPECT_THAT(ifnity->flight_plan().number_of_manœuvres(), Eq(16));
  std::vector<std::pair<DateTime, Speed>> manœuvre_ignition_tt_seconds_and_Δvs;
  for (int i = 0; i < ifnity->flight_plan().number_of_manœuvres(); ++i) {
    manœuvre_ignition_tt_seconds_and_Δvs.emplace_back(
        TTSecond(ifnity->flight_plan().GetManœuvre(i).initial_time()),
        ifnity->flight_plan().GetManœuvre(i).Δv().Norm());
  }
  // The flight plan only covers the inner solar system (this is probably
  // because of #3035).
  // It also differs from https://youtu.be/7BDxZV7UD9I?t=439.
  // TODO(egg): Compute the flybys and figure out what exactly is going on in
  // this flight plan.
  EXPECT_THAT(manœuvre_ignition_tt_seconds_and_Δvs,
              ElementsAre(Pair("1970-08-14T09:34:49"_DateTime,
                               3.80488671073918022e+03 * (Metre / Second)),
                          Pair("1970-08-15T13:59:24"_DateTime,
                               3.04867185471741759e-04 * (Metre / Second)),
                          Pair("1970-12-22T07:48:21"_DateTime,
                               1.58521291818444873e-03 * (Metre / Second)),
                          Pair("1971-01-08T17:36:55"_DateTime,
                               1.40000000034068623e-03 * (Metre / Second)),
                          Pair("1971-07-02T17:16:00"_DateTime,
                               1.00000000431022681e-04 * (Metre / Second)),
                          Pair("1971-09-06T03:27:33"_DateTime,
                               1.78421858738381537e-03 * (Metre / Second)),
                          Pair("1972-02-13T22:47:26"_DateTime,
                               7.72606625794511597e-04 * (Metre / Second)),
                          Pair("1972-03-25T16:30:19"_DateTime,
                               5.32846131747503372e-03 * (Metre / Second)),
                          Pair("1972-12-24T04:09:32"_DateTime,
                               3.45000000046532824e-03 * (Metre / Second)),
                          Pair("1973-06-04T01:59:07"_DateTime,
                               9.10695453328359134e-03 * (Metre / Second)),
                          Pair("1973-07-09T06:07:17"_DateTime,
                               4.49510921430966881e-01 * (Metre / Second)),
                          Pair("1973-09-10T03:59:44"_DateTime,
                               1.00000000431022681e-04 * (Metre / Second)),
                          Pair("1974-11-20T17:34:27"_DateTime,
                               5.10549409572428781e-01 * (Metre / Second)),
                          Pair("1975-10-07T01:29:45"_DateTime,
                               2.86686518692948443e-02 * (Metre / Second)),
                          Pair("1975-12-29T21:27:13"_DateTime,
                               1.00404183285598275e-03 * (Metre / Second)),
                          Pair("1977-07-28T22:47:53"_DateTime,
                               1.39666705839172456e-01 * (Metre / Second))));
}

}  // namespace ksp_plugin
}  // namespace principia
