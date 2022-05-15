
#include "physics/equipotential.hpp"

#include "base/not_null.hpp"
#include "geometry/frame.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/named_quantities.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "integrators/methods.hpp"
#include "integrators/embedded_explicit_runge_kutta_integrator.hpp"
#include "integrators/symmetric_linear_multistep_integrator.hpp"
#include "mathematica/mathematica.hpp"
#include "physics/body_centred_non_rotating_dynamic_frame.hpp"
#include "physics/ephemeris.hpp"
#include "physics/solar_system.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/solar_system_factory.hpp"

namespace principia {
namespace physics {

using base::make_not_null_unique;
using base::not_null;
using geometry::Arbitrary;
using geometry::Bivector;
using geometry::Frame;
using geometry::Inertial;
using geometry::Instant;
using geometry::Position;
using integrators::EmbeddedExplicitRungeKuttaIntegrator;
using integrators::SymmetricLinearMultistepIntegrator;
using integrators::methods::DormandPrince1986RK547FC;
using integrators::methods::QuinlanTremaine1990Order12;
using quantities::si::Hour;
using quantities::si::Metre;
using quantities::si::Milli;
using quantities::si::Minute;
using testing_utilities::SolarSystemFactory;

class EquipotentialTest : public ::testing::Test {
 protected:
  using Barycentric = Frame<enum class BarycentricTag, Inertial>;
  using World = Frame<enum class WorldTag, Arbitrary>;
  EquipotentialTest()
      : ephemeris_parameters_(
            SymmetricLinearMultistepIntegrator<QuinlanTremaine1990Order12,
                                               Position<Barycentric>>(),
            /*step=*/10 * Minute),
        solar_system_(make_not_null_unique<SolarSystem<Barycentric>>(
            SOLUTION_DIR / "astronomy" / "sol_gravity_model.proto.txt",
            SOLUTION_DIR / "astronomy" /
                "sol_initial_state_jd_2451545_000000000.proto.txt",
            /*ignore_frame=*/true)),
        ephemeris_(solar_system_->MakeEphemeris(
            /*accuracy_parameters=*/{/*fitting_tolerance=*/1 * Milli(Metre),
                                     /*geopotential_tolerance=*/0x1p-24},
            ephemeris_parameters_)),
        dynamic_frame_(BodyCentredNonRotatingDynamicFrame<Barycentric, World>(
            ephemeris_.get(),
            solar_system_->massive_body(
                *ephemeris_,
                SolarSystemFactory::name(SolarSystemFactory::Sun)))),
        equipotential_parameters_(
            EmbeddedExplicitRungeKuttaIntegrator<DormandPrince1986RK547FC,
                                                 Position<World>,
                                                 double>(),
            /*max_steps=*/1000,
            /*length_integration_tolerance=*/1 * Metre),
        equipotential_(equipotential_parameters_, &dynamic_frame_) {}

  Position<World> ComputePositionInWorld(SolarSystemFactory::Index const body) {
    auto const to_this_frame = dynamic_frame_.ToThisFrameAtTime(t0_);
    return to_this_frame.rigid_transformation()(
        solar_system_->trajectory(*ephemeris_, SolarSystemFactory::name(body))
            .EvaluatePosition(t0_));
  }

  Instant const t0_;
  Ephemeris<Barycentric>::FixedStepParameters const ephemeris_parameters_;
  not_null<std::unique_ptr<SolarSystem<Barycentric>>> const solar_system_;
  not_null<std::unique_ptr<Ephemeris<Barycentric>>> const ephemeris_;
  BodyCentredNonRotatingDynamicFrame<Barycentric, World> const dynamic_frame_;
  Equipotential<Barycentric, World>::AdaptiveParameters const
      equipotential_parameters_;
  Equipotential<Barycentric, World> equipotential_;
};

TEST_F(EquipotentialTest, Mathematica) {
  mathematica::Logger logger(TEMP_DIR / "equipotential.wl");
  Bivector<double, World> const plane({2, 3, -5});
  {
    LOG(ERROR)<<"MERCURY";
    Instant const t1 = t0_ + 24 * Hour;
    CHECK_OK(ephemeris_->Prolong(t1));
    auto const& [positions, βs] = equipotential_.ComputeLine(
        plane,
        ComputePositionInWorld(SolarSystemFactory::Mercury),
        t1);
    logger.Set(
        "equipotentialMercury", positions, mathematica::ExpressIn(Metre));
    logger.Set("betaMercury", βs);
  }
  {
    LOG(ERROR)<<"EARTH";
    Instant const t1 = t0_ + 24 * Hour;
    CHECK_OK(ephemeris_->Prolong(t1));
    auto const& [positions, βs] = equipotential_.ComputeLine(
        plane,
        ComputePositionInWorld(SolarSystemFactory::Earth),
        t1);
    logger.Set("equipotentialEarth", positions, mathematica::ExpressIn(Metre));
    logger.Set("betaEarth", βs);
  }
  {
    LOG(ERROR)<<"JUPITER CLOSE";
    Instant const t1 = t0_ + 24 * Hour;
    CHECK_OK(ephemeris_->Prolong(t1));
    auto const& [positions, βs] = equipotential_.ComputeLine(
        plane,
        ComputePositionInWorld(SolarSystemFactory::Jupiter),
        t1);
    logger.Set(
        "equipotentialJupiterClose", positions, mathematica::ExpressIn(Metre));
    logger.Set("betaJupiterClose", βs);
  }
  {
    LOG(ERROR)<<"JUPITER FAR";
    Instant const t1 = t0_ + 2400 * Hour;
    CHECK_OK(ephemeris_->Prolong(t1));
    auto const& [positions, βs] = equipotential_.ComputeLine(
        plane,
        ComputePositionInWorld(SolarSystemFactory::Jupiter),
        t1);
    logger.Set(
        "equipotentialJupiterFar", positions, mathematica::ExpressIn(Metre));
    logger.Set("betaJupiterFar", βs);
  }
}

}  // namespace physics
}  // namespace principia
