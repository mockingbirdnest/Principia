
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
#include "physics/solar_system.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/solar_system_factory.hpp"

namespace principia {
namespace physics {

using base::make_not_null_unique;
using base::not_null;
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
        equipotential_parameters_(
            EmbeddedExplicitRungeKuttaIntegrator<DormandPrince1986RK547FC,
                                                 Position<Barycentric>,
                                                 double>(),
            /*max_steps=*/1000,
            /*length_integration_tolerance=*/1 * Metre),
        equipotential_(equipotential_parameters_, *ephemeris_) {}

  Instant const t0_;
  Ephemeris<Barycentric>::FixedStepParameters const ephemeris_parameters_;
  not_null<std::unique_ptr<SolarSystem<Barycentric>>> const solar_system_;
  not_null<std::unique_ptr<Ephemeris<Barycentric>>> const ephemeris_;
  Equipotential<Barycentric>::AdaptiveParameters const
      equipotential_parameters_;
  Equipotential<Barycentric> equipotential_;
};

TEST_F(EquipotentialTest, Mathematica) {
  mathematica::Logger logger(TEMP_DIR / "equipotential.wl");
  Bivector<double, Barycentric> const plane({2, 3, -5});
  {
    LOG(ERROR)<<"MERCURY";
    Instant const t1 = t0_ + 24 * Hour;
    CHECK_OK(ephemeris_->Prolong(t1));
    auto const& [positions, βs] = equipotential_.ComputeLine(
        plane,
        solar_system_
            ->trajectory(*ephemeris_,
                         SolarSystemFactory::name(SolarSystemFactory::Mercury))
            .EvaluatePosition(t0_),
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
        solar_system_
            ->trajectory(*ephemeris_,
                         SolarSystemFactory::name(SolarSystemFactory::Earth))
            .EvaluatePosition(t0_),
        t1);
    logger.Set("equipotentialEarth", positions, mathematica::ExpressIn(Metre));
    logger.Set("betaEarth", βs);
  }
  {
    LOG(ERROR)<<"JUPITER";
    Instant const t1 = t0_ + 2400 * Hour;
    CHECK_OK(ephemeris_->Prolong(t1));
    auto const& [positions, βs] = equipotential_.ComputeLine(
        plane,
        solar_system_
            ->trajectory(*ephemeris_,
                         SolarSystemFactory::name(SolarSystemFactory::Jupiter))
            .EvaluatePosition(t0_),
        t1);
    logger.Set(
        "equipotentialJupiter", positions, mathematica::ExpressIn(Metre));
    logger.Set("betaJupiter", βs);
  }
}

}  // namespace physics
}  // namespace principia
