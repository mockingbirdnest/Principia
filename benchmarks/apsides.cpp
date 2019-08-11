
// .\Release\x64\benchmarks.exe --benchmark_repetitions=3 --benchmark_filter=Apsides --benchmark_min_time=30  // NOLINT(whitespace/line_length)

#include "astronomy/epoch.hpp"
#include "astronomy/frames.hpp"
#include "astronomy/standard_product_3.hpp"
#include "base/not_null.hpp"
#include "benchmark/benchmark.h"
#include "geometry/grassmann.hpp"
#include "geometry/named_quantities.hpp"
#include "integrators/embedded_explicit_runge_kutta_nyström_integrator.hpp"
#include "integrators/methods.hpp"
#include "integrators/symmetric_linear_multistep_integrator.hpp"
#include "physics/apsides.hpp"
#include "physics/body_centred_non_rotating_dynamic_frame.hpp"
#include "physics/body_surface_dynamic_frame.hpp"
#include "physics/ephemeris.hpp"
#include "physics/solar_system.hpp"
#include "quantities/astronomy.hpp"
#include "quantities/si.hpp"

namespace principia {

using astronomy::GCRS;
using astronomy::ICRS;
using astronomy::InfiniteFuture;
using astronomy::ITRS;
using astronomy::StandardProduct3;
using base::dynamic_cast_not_null;
using base::not_null;
using geometry::Position;
using geometry::Vector;
using integrators::EmbeddedExplicitRungeKuttaNyströmIntegrator;
using integrators::methods::DormandالمكاوىPrince1986RKN434FM;
using integrators::methods::QuinlanTremaine1990Order12;
using integrators::SymmetricLinearMultistepIntegrator;
using quantities::astronomy::JulianYear;
using quantities::si::Metre;
using quantities::si::Milli;
using quantities::si::Minute;
using quantities::si::Second;

namespace physics {

class ApsidesBenchmark : public benchmark::Fixture {
 protected:
  void SetUp(benchmark::State&) override {
    if (!setup_done_) {
      solar_system_2010_ = std::make_unique<SolarSystem<ICRS>>(
          SOLUTION_DIR / "astronomy" / "sol_gravity_model.proto.txt",
          SOLUTION_DIR / "astronomy" /
              "sol_initial_state_jd_2455200_500000000.proto.txt");
      ephemeris_ = solar_system_2010_->MakeEphemeris(
          /*accuracy_parameters=*/{/*fitting_tolerance=*/5 * Milli(Metre),
                                   /*geopotential_tolerance=*/0x1p-24},
          Ephemeris<ICRS>::FixedStepParameters(
              SymmetricLinearMultistepIntegrator<QuinlanTremaine1990Order12,
                                                 Position<ICRS>>(),
              /*step=*/10 * Minute));
      earth_ = dynamic_cast_not_null<OblateBody<ICRS> const*>(
          solar_system_2010_->massive_body(*ephemeris_, "Earth"));
      earth_trajectory_ = ephemeris_->trajectory(earth_);

      StandardProduct3 const ilrsa_lageos2_sp3(
          SOLUTION_DIR / "astronomy" / "standard_product_3" /
              "ilrsa.orb.lageos2.160319.v35.sp3",
          StandardProduct3::Dialect::ILRSA);
      StandardProduct3::SatelliteIdentifier const lageos2_id{
          StandardProduct3::SatelliteGroup::General, 52};

      auto const ilrsa_lageos2_trajectory_itrs =
          ilrsa_lageos2_sp3.orbit(lageos2_id).front();
      auto const begin = ilrsa_lageos2_trajectory_itrs->Begin();
      ephemeris_->Prolong(begin.time());

      BodySurfaceDynamicFrame<ICRS, ITRS> const itrs(ephemeris_.get(), earth_);
      ilrsa_lageos2_trajectory_icrs_.Append(
          begin.time(),
          itrs.FromThisFrameAtTime(begin.time())(begin.degrees_of_freedom()));
      ephemeris_->FlowWithAdaptiveStep(
          &ilrsa_lageos2_trajectory_icrs_,
          Ephemeris<ICRS>::NoIntrinsicAcceleration,
          begin.time() + 1 * JulianYear,
          Ephemeris<ICRS>::AdaptiveStepParameters(
              EmbeddedExplicitRungeKuttaNyströmIntegrator<
                  DormandالمكاوىPrince1986RKN434FM,
                  Position<ICRS>>(),
              std::numeric_limits<std::int64_t>::max(),
              /*length_integration_tolerance=*/1 * Milli(Metre),
              /*speed_integration_tolerance=*/1 * Milli(Metre) / Second),
          /*max_ephemeris_steps=*/std::numeric_limits<std::int64_t>::max(),
          /*last_point_only=*/false);

      BodyCentredNonRotatingDynamicFrame<ICRS, GCRS> const gcrs(
          ephemeris_.get(), earth_);
      for (auto it = ilrsa_lageos2_trajectory_icrs_.Begin();
           it != ilrsa_lageos2_trajectory_icrs_.End();
           ++it) {
        ilrsa_lageos2_trajectory_gcrs_.Append(
            it.time(),
            gcrs.ToThisFrameAtTime(it.time())(it.degrees_of_freedom()));
      }

      setup_done_ = true;
    }
  }

  bool setup_done_ = false;
  std::unique_ptr<SolarSystem<ICRS>> solar_system_2010_;
  std::unique_ptr<Ephemeris<ICRS>> ephemeris_;
  OblateBody<ICRS> const* earth_;
  ContinuousTrajectory<ICRS> const* earth_trajectory_;
  DiscreteTrajectory<ICRS> ilrsa_lageos2_trajectory_icrs_;
  DiscreteTrajectory<GCRS> ilrsa_lageos2_trajectory_gcrs_;
};

BENCHMARK_F(ApsidesBenchmark, ComputeApsides)(benchmark::State& state) {
  for (auto _ : state) {
    DiscreteTrajectory<ICRS> apoapsides;
    DiscreteTrajectory<ICRS> periapsides;
    ComputeApsides(*earth_trajectory_,
                   ilrsa_lageos2_trajectory_icrs_.Begin(),
                   ilrsa_lageos2_trajectory_icrs_.End(),
                   apoapsides,
                   periapsides);
    CHECK_EQ(2364, apoapsides.Size());
    CHECK_EQ(2365, periapsides.Size());
  }
}

BENCHMARK_F(ApsidesBenchmark, ComputeNodes)(benchmark::State& state) {
  for (auto _ : state) {
    DiscreteTrajectory<GCRS> ascending;
    DiscreteTrajectory<GCRS> descending;
    ComputeNodes(ilrsa_lageos2_trajectory_gcrs_.Begin(),
                 ilrsa_lageos2_trajectory_gcrs_.End(),
                 Vector<double, GCRS>({0, 0, 1}),
                 ascending,
                 descending);
    CHECK_EQ(2365, ascending.Size());
    CHECK_EQ(2365, descending.Size());
  }
}

}  // namespace physics
}  // namespace principia
