
#include "physics/equipotential.hpp"

#include <string>
#include <vector>

#include "absl/strings/str_cat.h"
#include "base/not_null.hpp"
#include "geometry/barycentre_calculator.hpp"
#include "geometry/frame.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/named_quantities.hpp"
#include "geometry/rotation.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "integrators/methods.hpp"
#include "integrators/embedded_explicit_runge_kutta_integrator.hpp"
#include "integrators/symmetric_linear_multistep_integrator.hpp"
#include "mathematica/mathematica.hpp"
#include "physics/body_centred_body_direction_dynamic_frame.hpp"
#include "physics/body_centred_non_rotating_dynamic_frame.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "physics/dynamic_frame.hpp"
#include "physics/ephemeris.hpp"
#include "physics/solar_system.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/solar_system_factory.hpp"

namespace principia {
namespace physics {

using base::make_not_null_unique;
using base::not_null;
using geometry::Arbitrary;
using geometry::Barycentre;
using geometry::Bivector;
using geometry::Frame;
using geometry::Inertial;
using geometry::Instant;
using geometry::Position;
using geometry::Rotation;
using geometry::Velocity;
using integrators::EmbeddedExplicitRungeKuttaIntegrator;
using integrators::SymmetricLinearMultistepIntegrator;
using integrators::methods::DormandPrince1986RK547FC;
using integrators::methods::QuinlanTremaine1990Order12;
using quantities::si::Day;
using quantities::si::Degree;
using quantities::si::Metre;
using quantities::si::Milli;
using quantities::si::Minute;
using quantities::si::Second;
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
        equipotential_parameters_(
            EmbeddedExplicitRungeKuttaIntegrator<
                DormandPrince1986RK547FC,
                Equipotential<Barycentric, World>::IndependentVariable,
                Position<World>,
                double>(),
            /*max_steps=*/1000,
            /*length_integration_tolerance=*/1 * Metre) {}

  Position<World> ComputePositionInWorld(
      Instant const& t,
      DynamicFrame<Barycentric, World> const& dynamic_frame,
      SolarSystemFactory::Index const body) {
    auto const to_this_frame = dynamic_frame.ToThisFrameAtTime(t);
    return to_this_frame.rigid_transformation()(
        solar_system_->trajectory(*ephemeris_, SolarSystemFactory::name(body))
            .EvaluatePosition(t));
  }

  // Logs to Mathematica the equipotential line for the given |body| in the
  // specified |dynamic_frame|.
  void LogEquipotentialLine(
      mathematica::Logger& logger,
      Bivector<double, World> const& plane,
      Instant const& t,
      DynamicFrame<Barycentric, World> const& dynamic_frame,
      SolarSystemFactory::Index const body,
      std::string_view const suffix = "") {
    Equipotential<Barycentric, World> const equipotential(
        equipotential_parameters_, &dynamic_frame);
    std::string const name = SolarSystemFactory::name(body);

    CHECK_OK(ephemeris_->Prolong(t));
    auto const& [positions, βs] = equipotential.ComputeLine(
        plane,
        t,
        ComputePositionInWorld(
            t0_, dynamic_frame, SolarSystemFactory::Mercury));
    logger.Set(absl::StrCat("equipotential", name, suffix),
               positions,
               mathematica::ExpressIn(Metre));
    logger.Set(absl::StrCat("beta", name, suffix), βs);
  }

  // Logs to Mathematica a family of equipotential lines determined by a
  // parameter.  There must exist an overload of |ComputeLine| with a
  // |LineParameter| as its third argument.
  template<typename LineParameter>
  void LogFamilyOfEquipotentialLines(
      mathematica::Logger& logger,
      DynamicFrame<Barycentric, World> const& dynamic_frame,
      int const number_of_days,
      std::string_view const suffix,
      std::function<std::vector<LineParameter>(
          Position<World> const& l4,
          Position<World> const& l5)> const& get_line_parameters) {
    Equipotential<Barycentric, World> const equipotential(
        equipotential_parameters_, &dynamic_frame);
    Bivector<double, World> const plane({0, 0, 1});

    std::vector<std::vector<std::vector<Position<World>>>> all_positions;
    std::vector<std::vector<std::vector<double>>> all_βs;
    for (int j = 0; j < number_of_days; ++j) {
      Instant const t = t0_ + j * Day;
      CHECK_OK(ephemeris_->Prolong(t));
      all_positions.emplace_back();
      all_βs.emplace_back();
      auto const earth_position =
          ComputePositionInWorld(t, dynamic_frame, SolarSystemFactory::Earth);
      auto const moon_position =
          ComputePositionInWorld(t, dynamic_frame, SolarSystemFactory::Moon);
      auto const moon_earth = earth_position - moon_position;

      Rotation<World, World> const rot_l4(-60 * Degree, plane);
      auto const moon_l4 = rot_l4(moon_earth);
      auto const l4 = moon_l4 + moon_position;
      Rotation<World, World> const rot_l5(60 * Degree, plane);
      auto const moon_l5 = rot_l5(moon_earth);
      auto const l5 = moon_l5 + moon_position;

      for (auto const& line_parameter : get_line_parameters(l4, l5)) {
        auto const& [positions, βs] =
            equipotential.ComputeLine(plane, t, line_parameter);
        all_positions.back().push_back(positions);
        all_βs.back().push_back(βs);
      }
    }
    logger.Set(absl::StrCat("equipotentialsEarthMoon", suffix),
               all_positions,
               mathematica::ExpressIn(Metre));
    logger.Set(absl::StrCat("betasEarthMoon", suffix), all_βs);
  }

  Instant const t0_;
  Ephemeris<Barycentric>::FixedStepParameters const ephemeris_parameters_;
  not_null<std::unique_ptr<SolarSystem<Barycentric>>> const solar_system_;
  not_null<std::unique_ptr<Ephemeris<Barycentric>>> const ephemeris_;
  Equipotential<Barycentric, World>::AdaptiveParameters const
      equipotential_parameters_;
};

#if !_DEBUG
TEST_F(EquipotentialTest, BodyCentredNonRotating) {
  mathematica::Logger logger(TEMP_DIR / "equipotential_bcnr.wl",
                             /*make_unique=*/false);
  auto const dynamic_frame(
      BodyCentredNonRotatingDynamicFrame<Barycentric, World>(
          ephemeris_.get(),
          solar_system_->massive_body(
              *ephemeris_, SolarSystemFactory::name(SolarSystemFactory::Sun))));
  Equipotential<Barycentric, World> const equipotential(
      equipotential_parameters_, &dynamic_frame);

  Bivector<double, World> const plane({2, 3, -5});

  LogEquipotentialLine(logger,
                       plane,
                       t0_ + 1 * Day,
                       dynamic_frame,
                       SolarSystemFactory::Mercury);
  LogEquipotentialLine(logger,
                       plane,
                       t0_ + 1 * Day,
                       dynamic_frame,
                       SolarSystemFactory::Earth);
  LogEquipotentialLine(logger,
                       plane,
                       t0_ + 1 * Day,
                       dynamic_frame,
                       SolarSystemFactory::Jupiter, "Close");
  LogEquipotentialLine(logger,
                       plane,
                       t0_ + 100 * Day,
                       dynamic_frame,
                       SolarSystemFactory::Jupiter, "Far");
}

TEST_F(EquipotentialTest, BodyCentredBodyDirection_EquidistantPoints) {
  mathematica::Logger logger(TEMP_DIR / "equipotential_bcbd.wl",
                             /*make_unique=*/true);
  auto const dynamic_frame(
      BodyCentredBodyDirectionDynamicFrame<Barycentric, World>(
          ephemeris_.get(),
          solar_system_->massive_body(
              *ephemeris_,
              SolarSystemFactory::name(SolarSystemFactory::Earth)),
          solar_system_->massive_body(
              *ephemeris_,
              SolarSystemFactory::name(SolarSystemFactory::Moon))));

  LogFamilyOfEquipotentialLines<Position<World>>(
      logger,
      dynamic_frame,
      /*number_of_days=*/30,
      /*suffix=*/"Positions",
      [](Position<World> const& l4, Position<World> const& l5) {
        std::vector<Position<World>> positions;
        for (int i = 0; i <= 10; ++i) {
          positions.push_back(Barycentre(
              std::pair{l4, l5}, std::pair{i / 10.0, (10.0 - i) / 10.0}));
        }
        return positions;
      });
}

TEST_F(EquipotentialTest, BodyCentredBodyDirection_EquidistantEnergies) {
  mathematica::Logger logger(TEMP_DIR / "equipotential_bcbd.wl",
                             /*make_unique=*/true);
  auto const dynamic_frame(
      BodyCentredBodyDirectionDynamicFrame<Barycentric, World>(
          ephemeris_.get(),
          solar_system_->massive_body(
              *ephemeris_,
              SolarSystemFactory::name(SolarSystemFactory::Earth)),
          solar_system_->massive_body(
              *ephemeris_,
              SolarSystemFactory::name(SolarSystemFactory::Moon))));

  LogFamilyOfEquipotentialLines<DegreesOfFreedom<World>>(
      logger,
      dynamic_frame,
      /*number_of_days=*/30,
      /*suffix=*/"Energies",
      [](Position<World> const& l4, Position<World> const& l5) {
        auto const midpoint = Barycentre(std::pair{l4, l5}, std::pair{1, 1});
        std::vector<DegreesOfFreedom<World>> degrees_of_freedom;
        for (int i = 0; i <= 10; ++i) {
          degrees_of_freedom.push_back(
              DegreesOfFreedom<World>(midpoint,
                                      Velocity<World>({i * 100 * Metre / Second,
                                                       0 * Metre / Second,
                                                       0 * Metre / Second})));
        }
        for (int i = 0; i <= 10; ++i) {
          degrees_of_freedom.push_back(DegreesOfFreedom<World>(
              midpoint,
              Velocity<World>({(1100 + i * 10) * Metre / Second,
                               0 * Metre / Second,
                               0 * Metre / Second})));
        }
        return degrees_of_freedom;
      });
}

#endif

}  // namespace physics
}  // namespace principia
