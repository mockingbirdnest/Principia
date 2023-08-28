// .\Release\x64\benchmarks.exe --benchmark_filter=LagrangeEquipotential

#include <vector>

#include "base/not_null.hpp"
#include "base/status_utilities.hpp"  // 🧙 For CHECK_OK.
#include "benchmark/benchmark.h"
#include "geometry/frame.hpp"
#include "geometry/instant.hpp"
#include "integrators/methods.hpp"
#include "integrators/symmetric_linear_multistep_integrator.hpp"
#include "physics/ephemeris.hpp"
#include "physics/lagrange_equipotentials.hpp"
#include "physics/massive_body.hpp"
#include "physics/solar_system.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/solar_system_factory.hpp"

namespace principia {
namespace physics {

using namespace principia::base::_not_null;
using namespace principia::geometry::_frame;
using namespace principia::geometry::_instant;
using namespace principia::integrators::_methods;
using namespace principia::integrators::_symmetric_linear_multistep_integrator;
using namespace principia::physics::_ephemeris;
using namespace principia::physics::_lagrange_equipotentials;
using namespace principia::physics::_massive_body;
using namespace principia::physics::_solar_system;
using namespace principia::quantities::_si;
using namespace principia::testing_utilities::_solar_system_factory;

using Barycentric = Frame<struct BarycentricTag, Inertial>;
using World = Frame<struct WorldTag, Arbitrary>;

class LagrangeEquipotentialsBenchmark : public benchmark::Fixture {
 protected:
  static void SetUpFixture() {
    Ephemeris<Barycentric>::FixedStepParameters const ephemeris_parameters(
        SymmetricLinearMultistepIntegrator<
            QuinlanTremaine1990Order12,
            Ephemeris<Barycentric>::NewtonianMotionEquation>(),
        /*step=*/10 * Minute);
    solar_system_ =
        make_not_null_unique<SolarSystem<Barycentric>>(
            SOLUTION_DIR / "astronomy" / "sol_gravity_model.proto.txt",
            SOLUTION_DIR / "astronomy" /
                "sol_initial_state_jd_2451545_000000000.proto.txt",
            /*ignore_frame=*/true)
            .release();
    ephemeris_ =
        solar_system_
            ->MakeEphemeris(
                /*accuracy_parameters=*/{/*fitting_tolerance=*/1 * Milli(Metre),
                                         /*geopotential_tolerance=*/0x1p-24},
                ephemeris_parameters)
            .release();
    CHECK_OK(ephemeris_->Prolong(t0_));
  }

  void SetUp(benchmark::State&) override {
    static int const set_up_fixture = []() {
      SetUpFixture();
      return 0;
    }();
  }

  static constexpr Instant t0_{};
  static SolarSystem<Barycentric>* solar_system_;
  static Ephemeris<Barycentric>* ephemeris_;
};

SolarSystem<Barycentric>* LagrangeEquipotentialsBenchmark::solar_system_;
Ephemeris<Barycentric>* LagrangeEquipotentialsBenchmark::ephemeris_;

BENCHMARK_F(LagrangeEquipotentialsBenchmark, EarthMoon)(
    benchmark::State& state) {
  auto const earth = solar_system_->massive_body(
      *ephemeris_, SolarSystemFactory::name(SolarSystemFactory::Earth));
  auto const moon = solar_system_->massive_body(
      *ephemeris_, SolarSystemFactory::name(SolarSystemFactory::Moon));

  for (auto _ : state) {
    auto const equipotentials =
        LagrangeEquipotentials<Barycentric, World>(ephemeris_)
            .ComputeLines({.primaries = {earth},
                           .secondaries = {moon},
                           .time = t0_});
    CHECK_OK(equipotentials.status());
  }
}

BENCHMARK_F(LagrangeEquipotentialsBenchmark, SunNeptune)(
    benchmark::State& state) {
  LagrangeEquipotentials<Barycentric, World>::Parameters parameters;
  for (int i = SolarSystemFactory::Sun; i <= SolarSystemFactory::LastBody;
       ++i) {
    switch (i) {
      case SolarSystemFactory::Pluto:
      case SolarSystemFactory::Charon:
      case SolarSystemFactory::Eris:
        continue;
      case SolarSystemFactory::Neptune:
      case SolarSystemFactory::Triton:
        parameters.secondaries.push_back(solar_system_->massive_body(
          *ephemeris_, SolarSystemFactory::name(i)));
        break;
      default:
        parameters.primaries.push_back(solar_system_->massive_body(
          *ephemeris_, SolarSystemFactory::name(i)));
    }
  }
  std::vector<not_null<MassiveBody const*>> primaries;

  for (auto _ : state) {
    parameters.time = t0_;
    auto const equipotentials =
        LagrangeEquipotentials<Barycentric, World>(ephemeris_)
            .ComputeLines(parameters);
    CHECK_OK(equipotentials.status());
  }
}

}  // namespace physics
}  // namespace principia
