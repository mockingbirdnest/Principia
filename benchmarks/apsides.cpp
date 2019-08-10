
#include "astronomy/frames.hpp"
#include "astronomy/standard_product_3.hpp"
#include "base/not_null.hpp"
#include "benchmark/benchmark.h"
#include "geometry/named_quantities.hpp"
#include "integrators/methods.hpp"
#include "integrators/symmetric_linear_multistep_integrator.hpp"
#include "physics/apsides.hpp"
#include "physics/body_surface_dynamic_frame.hpp"
#include "physics/ephemeris.hpp"
#include "physics/solar_system.hpp"
#include "quantities/si.hpp"

namespace principia {

using astronomy::ICRS;
using astronomy::ITRS;
using astronomy::StandardProduct3;
using base::dynamic_cast_not_null;
using base::not_null;
using geometry::Position;
using integrators::methods::QuinlanTremaine1990Order12;
using integrators::SymmetricLinearMultistepIntegrator;
using quantities::si::Metre;
using quantities::si::Milli;
using quantities::si::Minute;

namespace physics {

class ApsidesBenchmark : public benchmark::Fixture {
 protected:
  void SetUp(const benchmark::State&) override {
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

    BodySurfaceDynamicFrame<ICRS, ITRS> const itrs(ephemeris_.get(), earth_);
    auto const ilrsa_lageos2_trajectory_itrs =
        ilrsa_lageos2_sp3.orbit(lageos2_id).front();
    ephemeris_->Prolong(ilrsa_lageos2_trajectory_itrs->last().time());
    for (auto it = ilrsa_lageos2_trajectory_itrs->Begin();
         it != ilrsa_lageos2_trajectory_itrs->End();
         ++it) {
      ilrsa_lageos2_trajectory_.Append(
          it.time(),
          itrs.FromThisFrameAtTime(it.time())(it.degrees_of_freedom()));
    }
  }

  std::unique_ptr<SolarSystem<ICRS>> solar_system_2010_;
  std::unique_ptr<Ephemeris<ICRS>> ephemeris_;
  OblateBody<ICRS> const* earth_;
  ContinuousTrajectory<ICRS> const* earth_trajectory_;
  DiscreteTrajectory<ICRS> ilrsa_lageos2_trajectory_;
};

BENCHMARK_F(ApsidesBenchmark, ComputeApsides)(benchmark::State& state) {
  for (auto _ : state) {
    DiscreteTrajectory<ICRS> apoapsides;
    DiscreteTrajectory<ICRS> periapsides;
    ComputeApsides(*earth_trajectory_,
                   ilrsa_lageos2_trajectory_.Begin(),
                   ilrsa_lageos2_trajectory_.End(),
                   apoapsides,
                   periapsides);
  }
}

}  // namespace physics
}  // namespace principia
