// .\Release\x64\benchmarks.exe --benchmark_repetitions=3 --benchmark_filter=Apsides --benchmark_min_time=30  // NOLINT(whitespace/line_length)

#include <limits>

#include "astronomy/epoch.hpp"
#include "astronomy/frames.hpp"
#include "astronomy/standard_product_3.hpp"
#include "base/not_null.hpp"
#include "base/status_utilities.hpp"
#include "benchmark/benchmark.h"
#include "geometry/grassmann.hpp"
#include "integrators/embedded_explicit_runge_kutta_nyström_integrator.hpp"
#include "integrators/methods.hpp"
#include "integrators/symmetric_linear_multistep_integrator.hpp"
#include "physics/apsides.hpp"
#include "physics/body_centred_non_rotating_reference_frame.hpp"
#include "physics/body_surface_reference_frame.hpp"
#include "physics/discrete_trajectory.hpp"
#include "physics/ephemeris.hpp"
#include "physics/solar_system.hpp"
#include "quantities/astronomy.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace physics {

using namespace principia::astronomy::_frames;
using namespace principia::astronomy::_standard_product_3;
using namespace principia::base::_not_null;
using namespace principia::geometry::_grassmann;
using namespace principia::integrators::_embedded_explicit_runge_kutta_nyström_integrator;  // NOLINT
using namespace principia::integrators::_methods;
using namespace principia::integrators::_symmetric_linear_multistep_integrator;
using namespace principia::quantities::_astronomy;
using namespace principia::quantities::_si;

class ApsidesBenchmark : public benchmark::Fixture {
 protected:
  // Benchmark doesn't have that capability, so we have to do it ourselves.
  // This function takes about 40 s to run.
  static void SetUpFixture() {
    solar_system_2010_ = std::make_unique<SolarSystem<ICRS>>(
        SOLUTION_DIR / "astronomy" / "sol_gravity_model.proto.txt",
        SOLUTION_DIR / "astronomy" /
            "sol_initial_state_jd_2455200_500000000.proto.txt").release();
    ephemeris_ =
        solar_system_2010_
            ->MakeEphemeris(
                /*accuracy_parameters=*/{/*fitting_tolerance=*/5 * Milli(Metre),
                                         /*geopotential_tolerance=*/0x1p-24},
                Ephemeris<ICRS>::FixedStepParameters(
                    SymmetricLinearMultistepIntegrator<
                        QuinlanTremaine1990Order12,
                        Ephemeris<ICRS>::NewtonianMotionEquation>(),
                    /*step=*/10 * Minute))
            .release();
    earth_ = dynamic_cast_not_null<OblateBody<ICRS> const*>(
        solar_system_2010_->massive_body(*ephemeris_, "Earth"));
    earth_trajectory_ = ephemeris_->trajectory(earth_);

    StandardProduct3 const ilrsa_lageos2_sp3(
        SOLUTION_DIR / "astronomy" / "standard_product_3" /
            "ilrsa.orb.lageos2.160319.v35.sp3",
        StandardProduct3::Dialect::ILRSA);
    CHECK(ilrsa_lageos2_sp3.file_has_velocities());
    StandardProduct3::SatelliteIdentifier const lageos2_id{
        StandardProduct3::SatelliteGroup::General, 52};

    auto const ilrsa_lageos2_trajectory_itrs =
        ilrsa_lageos2_sp3.orbit(lageos2_id).front();
    auto const begin = ilrsa_lageos2_trajectory_itrs->begin();
    CHECK_OK(ephemeris_->Prolong(begin->time));

    BodySurfaceReferenceFrame<ICRS, ITRS> const itrs(ephemeris_, earth_);
    ilrsa_lageos2_trajectory_icrs_ = new DiscreteTrajectory<ICRS>;
    CHECK_OK(ilrsa_lageos2_trajectory_icrs_->Append(
        begin->time,
        itrs.FromThisFrameAtTime(begin->time)(begin->degrees_of_freedom)));
    CHECK_OK(ephemeris_->FlowWithAdaptiveStep(
        ilrsa_lageos2_trajectory_icrs_,
        Ephemeris<ICRS>::NoIntrinsicAcceleration,
        begin->time + 1 * JulianYear,
        Ephemeris<ICRS>::AdaptiveStepParameters(
            EmbeddedExplicitRungeKuttaNyströmIntegrator<
                DormandالمكاوىPrince1986RKN434FM,
                Ephemeris<ICRS>::NewtonianMotionEquation>(),
            std::numeric_limits<std::int64_t>::max(),
            /*length_integration_tolerance=*/1 * Milli(Metre),
            /*speed_integration_tolerance=*/1 * Milli(Metre) / Second),
        /*max_ephemeris_steps=*/std::numeric_limits<std::int64_t>::max()));

    BodyCentredNonRotatingReferenceFrame<ICRS, GCRS> const gcrs(ephemeris_,
                                                              earth_);
    ilrsa_lageos2_trajectory_gcrs_ = new DiscreteTrajectory<GCRS>;
    for (auto const& [time, degrees_of_freedom] :
         *ilrsa_lageos2_trajectory_icrs_) {
      CHECK_OK(ilrsa_lageos2_trajectory_gcrs_->Append(
          time, gcrs.ToThisFrameAtTime(time)(degrees_of_freedom)));
    }
  }

  void SetUp(benchmark::State&) override {
    static int const set_up_fixture = []() {
      SetUpFixture();
      return 0;
    }();
  }

  static SolarSystem<ICRS>* solar_system_2010_;
  static Ephemeris<ICRS>* ephemeris_;
  static OblateBody<ICRS> const* earth_;
  static ContinuousTrajectory<ICRS> const* earth_trajectory_;
  static DiscreteTrajectory<ICRS>* ilrsa_lageos2_trajectory_icrs_;
  static DiscreteTrajectory<GCRS>* ilrsa_lageos2_trajectory_gcrs_;
};

SolarSystem<ICRS>* ApsidesBenchmark::solar_system_2010_ = nullptr;
Ephemeris<ICRS>* ApsidesBenchmark::ephemeris_ = nullptr;
OblateBody<ICRS> const* ApsidesBenchmark::earth_ = nullptr;
ContinuousTrajectory<ICRS> const* ApsidesBenchmark::earth_trajectory_ = nullptr;
DiscreteTrajectory<ICRS>* ApsidesBenchmark::ilrsa_lageos2_trajectory_icrs_ =
    nullptr;
DiscreteTrajectory<GCRS>* ApsidesBenchmark::ilrsa_lageos2_trajectory_gcrs_ =
    nullptr;

BENCHMARK_F(ApsidesBenchmark, ComputeApsides)(benchmark::State& state) {
  for (auto _ : state) {
    DiscreteTrajectory<ICRS> apoapsides;
    DiscreteTrajectory<ICRS> periapsides;
    ComputeApsides(*earth_trajectory_,
                   *ilrsa_lageos2_trajectory_icrs_,
                   ilrsa_lageos2_trajectory_icrs_->begin(),
                   ilrsa_lageos2_trajectory_icrs_->end(),
                   /*max_points=*/std::numeric_limits<int>::max(),
                   apoapsides,
                   periapsides);
    CHECK_EQ(2364, apoapsides.size());
    CHECK_EQ(2365, periapsides.size());
  }
}

BENCHMARK_F(ApsidesBenchmark, ComputeNodes)(benchmark::State& state) {
  for (auto _ : state) {
    DiscreteTrajectory<GCRS> ascending;
    DiscreteTrajectory<GCRS> descending;
    CHECK_OK(ComputeNodes(*ilrsa_lageos2_trajectory_gcrs_,
                          ilrsa_lageos2_trajectory_gcrs_->begin(),
                          ilrsa_lageos2_trajectory_gcrs_->end(),
                          Vector<double, GCRS>({0, 0, 1}),
                          /*max_points=*/std::numeric_limits<int>::max(),
                          ascending,
                          descending));
    CHECK_EQ(2365, ascending.size());
    CHECK_EQ(2365, descending.size());
  }
}

}  // namespace physics
}  // namespace principia
