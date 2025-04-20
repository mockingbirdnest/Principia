#include "physics/equipotential.hpp"

#include <array>
#include <functional>
#include <memory>
#include <string>
#include <vector>

#include "absl/strings/str_cat.h"
#include "astronomy/stabilize_ksp.hpp"
#include "base/not_null.hpp"
#include "geometry/barycentre_calculator.hpp"
#include "geometry/frame.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/instant.hpp"
#include "geometry/plane.hpp"
#include "geometry/rotation.hpp"
#include "geometry/space.hpp"
#include "gtest/gtest.h"
#include "integrators/embedded_explicit_runge_kutta_integrator.hpp"
#include "integrators/methods.hpp"
#include "integrators/symmetric_linear_multistep_integrator.hpp"
#include "integrators/symplectic_runge_kutta_nyström_integrator.hpp"
#include "mathematica/logger.hpp"
#include "mathematica/mathematica.hpp"
#include "physics/body_centred_body_direction_reference_frame.hpp"
#include "physics/body_centred_non_rotating_reference_frame.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "physics/ephemeris.hpp"
#include "physics/reference_frame.hpp"
#include "physics/rotating_pulsating_reference_frame.hpp"
#include "physics/solar_system.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/solar_system_factory.hpp"

namespace principia {
namespace physics {

using namespace principia::astronomy::_stabilize_ksp;
using namespace principia::base::_not_null;
using namespace principia::geometry::_barycentre_calculator;
using namespace principia::geometry::_frame;
using namespace principia::geometry::_grassmann;
using namespace principia::geometry::_instant;
using namespace principia::geometry::_plane;
using namespace principia::geometry::_rotation;
using namespace principia::geometry::_space;
using namespace principia::integrators::_embedded_explicit_runge_kutta_integrator;  // NOLINT
using namespace principia::integrators::_methods;
using namespace principia::integrators::_symmetric_linear_multistep_integrator;
using namespace principia::integrators::_symplectic_runge_kutta_nyström_integrator;  // NOLINT
using namespace principia::mathematica::_logger;
using namespace principia::mathematica::_mathematica;
using namespace principia::physics::_body_centred_body_direction_reference_frame;  // NOLINT
using namespace principia::physics::_body_centred_non_rotating_reference_frame;
using namespace principia::physics::_degrees_of_freedom;
using namespace principia::physics::_ephemeris;
using namespace principia::physics::_equipotential;
using namespace principia::physics::_reference_frame;
using namespace principia::physics::_rotating_pulsating_reference_frame;
using namespace principia::physics::_solar_system;
using namespace principia::quantities::_named_quantities;
using namespace principia::quantities::_quantities;
using namespace principia::quantities::_si;
using namespace principia::testing_utilities::_solar_system_factory;

class EquipotentialTest : public ::testing::Test {
 protected:
  using Barycentric = Frame<struct BarycentricTag, Inertial>;
  using World = Frame<struct WorldTag, Arbitrary>;

  EquipotentialTest()
      : ephemeris_parameters_(
            SymmetricLinearMultistepIntegrator<
                QuinlanTremaine1990Order12,
                Ephemeris<Barycentric>::NewtonianMotionEquation>(),
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
                Equipotential<Barycentric, World>::ODE>(),
            /*max_steps=*/1000,
            /*length_integration_tolerance=*/1 * Metre) {}

  Position<World> ComputePositionInWorld(
      Instant const& t,
      ReferenceFrame<Barycentric, World> const& reference_frame,
      SolarSystemFactory::Index const body) {
    auto const to_this_frame = reference_frame.ToThisFrameAtTimeSimilarly(t);
    return to_this_frame.similarity()(
        solar_system_->trajectory(*ephemeris_, SolarSystemFactory::name(body))
            .EvaluatePosition(t));
  }

  std::array<Position<World>, 2> ComputeLagrangePoints(
      SolarSystemFactory::Index const body1,
      SolarSystemFactory::Index const body2,
      Instant const& t,
      ReferenceFrame<Barycentric, World> const& reference_frame,
      Plane<World> const& plane) {
    auto const body1_position =
        ComputePositionInWorld(t, reference_frame, body1);
    auto const body2_position =
        ComputePositionInWorld(t, reference_frame, body2);
    auto const body2_body1 = body1_position - body2_position;

    auto const& binormal = plane.UnitBinormals().front();
    Rotation<World, World> const rot_l4(-60 * Degree, binormal);
    auto const body2_l4 = rot_l4(body2_body1);
    auto const l4 = body2_l4 + body2_position;
    Rotation<World, World> const rot_l5(60 * Degree, binormal);
    auto const body2_l5 = rot_l5(body2_body1);
    auto const l5 = body2_l5 + body2_position;

    return {l4, l5};
  }

  // Logs to Mathematica the equipotential line for the given `body` in the
  // specified `reference_frame`.
  void LogEquipotentialLine(
      Logger& logger,
      Plane<World> const& plane,
      Instant const& t,
      ReferenceFrame<Barycentric, World> const& reference_frame,
      SolarSystemFactory::Index const body,
      std::string_view const suffix = "") {
    Equipotential<Barycentric, World> const equipotential(
        equipotential_parameters_,
        &reference_frame,
        /*characteristic_length=*/1 * Metre);
    std::string const name = SolarSystemFactory::name(body);

    CHECK_OK(ephemeris_->Prolong(t));
    auto const line =
        equipotential.ComputeLine(
            plane, t, ComputePositionInWorld(t0_, reference_frame, body));
    std::vector<Position<World>> positions;
    for (auto const& [s, dof] : line) {
      positions.push_back(dof.position());
    }
    logger.Set(absl::StrCat("equipotential", name, suffix),
               positions,
               ExpressIn(Metre));
  }

  // Logs to Mathematica a family of equipotential lines determined by a
  // parameter.  There must exist an overload of `ComputeLine` with a
  // `LineParameter` as its third argument.
  template<typename LineParameter>
  void LogFamilyOfEquipotentialLines(
      Logger& logger,
      ReferenceFrame<Barycentric, World> const& reference_frame,
      int const number_of_days,
      std::string_view const suffix,
      std::function<std::vector<LineParameter>(
          Position<World> const& l4,
          Position<World> const& l5)> const& get_line_parameters) {
    Equipotential<Barycentric, World> const equipotential(
        equipotential_parameters_,
        &reference_frame,
        /*characteristic_length=*/1 * Metre);
    auto const plane =
        Plane<World>::OrthogonalTo(Vector<double, World>({0, 0, 1}));

    std::vector<std::vector<std::vector<Position<World>>>> all_positions;
    for (int j = 0; j < number_of_days; ++j) {
      Instant const t = t0_ + j * Day;
      CHECK_OK(ephemeris_->Prolong(t));
      all_positions.emplace_back();

      auto const& [l4, l5] = ComputeLagrangePoints(SolarSystemFactory::Earth,
                                                   SolarSystemFactory::Moon,
                                                   t,
                                                   reference_frame,
                                                   plane);

      for (auto const& line_parameter : get_line_parameters(l4, l5)) {
        auto const line =
            equipotential.ComputeLine(plane, t, line_parameter);
        all_positions.back().emplace_back();
        for (auto const& [s, dof] : line) {
          all_positions.back().back().push_back(dof.position());
        }
      }
    }
    logger.Set(absl::StrCat("equipotentialsEarthMoon", suffix),
               all_positions,
               ExpressIn(Metre));
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
  Logger logger(TEMP_DIR / "equipotential_bcnr.wl",
                /*make_unique=*/false);
  auto const reference_frame(
      BodyCentredNonRotatingReferenceFrame<Barycentric, World>(
          ephemeris_.get(),
          solar_system_->massive_body(
              *ephemeris_, SolarSystemFactory::name(SolarSystemFactory::Sun))));
  Equipotential<Barycentric, World> const equipotential(
      equipotential_parameters_,
      &reference_frame,
      /*characteristic_length=*/1 * Metre);

  auto const plane =
      Plane<World>::OrthogonalTo(Vector<double, World>({2, 3, -5}));

  LogEquipotentialLine(logger,
                       plane,
                       t0_ + 1 * Day,
                       reference_frame,
                       SolarSystemFactory::Mercury);
  LogEquipotentialLine(logger,
                       plane,
                       t0_ + 1 * Day,
                       reference_frame,
                       SolarSystemFactory::Earth);
  LogEquipotentialLine(logger,
                       plane,
                       t0_ + 1 * Day,
                       reference_frame,
                       SolarSystemFactory::Jupiter, "Close");
  LogEquipotentialLine(logger,
                       plane,
                       t0_ + 100 * Day,
                       reference_frame,
                       SolarSystemFactory::Jupiter, "Far");
}

TEST_F(EquipotentialTest, BodyCentredBodyDirection_EquidistantPoints) {
  Logger logger(TEMP_DIR / "equipotential_bcbd_distances.wl",
                /*make_unique=*/false);
  auto const reference_frame(
      BodyCentredBodyDirectionReferenceFrame<Barycentric, World>(
          ephemeris_.get(),
          solar_system_->massive_body(
              *ephemeris_,
              SolarSystemFactory::name(SolarSystemFactory::Earth)),
          solar_system_->massive_body(
              *ephemeris_,
              SolarSystemFactory::name(SolarSystemFactory::Moon))));

  LogFamilyOfEquipotentialLines<Position<World>>(
      logger,
      reference_frame,
      /*number_of_days=*/30,
      /*suffix=*/"Distances",
      [](Position<World> const& l4, Position<World> const& l5) {
        std::vector<Position<World>> positions;
        for (int i = 0; i <= 10; ++i) {
          positions.push_back(
              Barycentre({l4, l5}, {i / 10.0, (10.0 - i) / 10.0}));
        }
        return positions;
      });
}
#endif

TEST_F(EquipotentialTest, KSPL₄L₅) {
  using KSP = Frame<struct KSPTag, Inertial>;
  using RotatingPulsating = Frame<struct RotatingPulsatingTag, Arbitrary>;

  constexpr Length box_side = 3 * Metre;
  constexpr Length characteristic_length = 1 * Nano(Metre);
  constexpr Length length_integration_tolerance = 10 * Femto(Metre);
  constexpr std::int64_t max_steps = 10'000;

  SolarSystem<KSP> solar_system(
      SOLUTION_DIR / "astronomy" / "kerbol_gravity_model.proto.txt",
      SOLUTION_DIR / "astronomy" / "kerbol_initial_state_0_0.proto.txt");
  StabilizeKSP(solar_system);

  auto const ephemeris = solar_system.MakeEphemeris(
      /*accuracy_parameters=*/{/*fitting_tolerance=*/1 * Milli(Metre),
                               /*geopotential_tolerance*/ 0x1.0p-24},
      Ephemeris<KSP>::FixedStepParameters(
          SymplecticRungeKuttaNyströmIntegrator<
              BlanesMoan2002SRKN14A,
              Ephemeris<KSP>::NewtonianMotionEquation>(),
          /*step=*/35 * Minute));

  auto const eve = solar_system.massive_body(*ephemeris, "Eve");
  auto const gilly = solar_system.massive_body(*ephemeris, "Gilly");
  auto const kerbin = solar_system.massive_body(*ephemeris, "Kerbin");
  auto const minmus = solar_system.massive_body(*ephemeris, "Minmus");
  auto const moho = solar_system.massive_body(*ephemeris, "Moho");
  auto const mun = solar_system.massive_body(*ephemeris, "Mun");
  auto const sun = solar_system.massive_body(*ephemeris, "Sun");
  auto const reference_frame(
      RotatingPulsatingReferenceFrame<KSP, RotatingPulsating>(
          ephemeris.get(), {sun, moho, eve, gilly}, {kerbin, mun, minmus}));

  Equipotential<KSP, RotatingPulsating> const equipotential(
      {EmbeddedExplicitRungeKuttaIntegrator<
           DormandPrince1986RK547FC,
           typename Equipotential<KSP, RotatingPulsating>::ODE>(),
       max_steps,
       length_integration_tolerance},
      &reference_frame,
      characteristic_length);

  auto const plane = Plane<RotatingPulsating>::OrthogonalTo(
      Vector<double, RotatingPulsating>({0, 0, 1}));
  Instant const t = Instant{} + 313.76000000000460 * Second;
  std::vector<Position<RotatingPulsating>> const peaks = {
      RotatingPulsating::origin +
          Displacement<RotatingPulsating>({-0.64802406026529646 * Metre,
                                           -0.76141434216037418 * Metre,
                                           0 * Metre}),
      RotatingPulsating::origin +
          Displacement<RotatingPulsating>({0.69811857380647024 * Metre,
                                           -0.71592199338059159 * Metre,
                                           0 * Metre}),
      RotatingPulsating::origin +
          Displacement<RotatingPulsating>({0.50997091222800128 * Metre,
                                           0.85996006756134868 * Metre,
                                           0 * Metre}),
      RotatingPulsating::origin +
          Displacement<RotatingPulsating>({-0.072379219947440787 * Metre,
                                           -0.99744151614018683 * Metre,
                                           0 * Metre}),
      RotatingPulsating::origin +
          Displacement<RotatingPulsating>({0.99692132254698529 * Metre,
                                           0.077369574947032782 * Metre,
                                           0 * Metre}),
      RotatingPulsating::origin +
          Displacement<RotatingPulsating>({-0.68469290995389276 * Metre,
                                           0.72835422306333097 * Metre,
                                           0 * Metre}),
      RotatingPulsating::origin +
          Displacement<RotatingPulsating>({-0.99958277597033018 * Metre,
                                           -0.036821342463918334 * Metre,
                                           0 * Metre}),
  };
  std::vector<Equipotential<KSP, RotatingPulsating>::Well> const wells = {
      {.position =
           RotatingPulsating::origin +
           Displacement<RotatingPulsating>({0.99999693050679217 * Metre,
                                            2.7392069189813769e-19 * Metre,
                                            -5.2603338992119579e-22 * Metre}),
       .radius = 0.0034803037593885916 * Metre},
      {.position =
           RotatingPulsating::origin +
           Displacement<RotatingPulsating>({-3.0694932080578120e-06 * Metre,
                                            -9.5558346840227855e-21 * Metre,
                                            -5.2603338992261777e-22 * Metre}),
       .radius = 0.73230345545398934 * Metre}};
  auto const towards_infinity = [](Position<RotatingPulsating> q) {
    return RotatingPulsating::origin +
           Normalize(q - RotatingPulsating::origin) * box_side;
  };
  SpecificEnergy const energy =
      -6.9913162881136685e-13 * si::Unit<SpecificEnergy>;

  CHECK_OK(ephemeris->Prolong(t));
  auto const lines = equipotential.ComputeLines(
      plane, t, peaks, wells, towards_infinity, energy);
}

}  // namespace physics
}  // namespace principia
