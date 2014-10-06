#include "physics/n_body_system.hpp"

#include <map>
#include <memory>
#include <string>
#include <vector>

#include "geometry/grassmann.hpp"
#include "geometry/named_quantities.hpp"
#include "geometry/point.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "physics/body.hpp"
#include "physics/trajectory.hpp"
#include "quantities/constants.hpp"
#include "quantities/numbers.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/almost_equals.hpp"
#include "testing_utilities/death_message.hpp"
#include "testing_utilities/numerics.hpp"
#include "testing_utilities/solar_system.hpp"

using principia::constants::GravitationalConstant;
using principia::geometry::Barycentre;
using principia::geometry::Instant;
using principia::geometry::Point;
using principia::geometry::Vector;
using principia::quantities::Angle;
using principia::quantities::ArcTan;
using principia::quantities::Area;
using principia::quantities::Pow;
using principia::quantities::SIUnit;
using principia::testing_utilities::AlmostEquals;
using principia::testing_utilities::DeathMessage;
using principia::testing_utilities::ICRFJ2000Ecliptic;
using principia::testing_utilities::kSolarSystemBarycentre;
using principia::testing_utilities::RelativeError;
using principia::testing_utilities::SolarSystem;
using principia::si::Degree;
using principia::si::Minute;
using testing::Eq;
using testing::Gt;
using testing::Lt;

namespace principia {
namespace physics {

class NBodySystemTest : public testing::Test {
 protected:
  struct EarthMoonOrbitPlane;

  void SetUp() override {
    integrator_.Initialize(integrator_.Order5Optimal());

    // The Earth-Moon system, roughly, with a circular orbit with velocities
    // in the center-of-mass frame.
    body1_ = std::make_unique<Body>(6E24 * SIUnit<Mass>());
    body2_ = std::make_unique<Body>(7E22 * SIUnit<Mass>());

    // A massless probe.
    body3_ = std::make_unique<Body>(0 * SIUnit<Mass>());

    trajectory1_ =
        std::make_unique<Trajectory<EarthMoonOrbitPlane>>(*body1_);
    trajectory2_ =
        std::make_unique<Trajectory<EarthMoonOrbitPlane>>(*body2_);
    trajectory3_ =
        std::make_unique<Trajectory<EarthMoonOrbitPlane>>(*body3_);
    Position<EarthMoonOrbitPlane> const q1(
        Vector<Length, EarthMoonOrbitPlane>({0 * SIUnit<Length>(),
                                             0 * SIUnit<Length>(),
                                             0 * SIUnit<Length>()}));
    Position<EarthMoonOrbitPlane> const q2(
        Vector<Length, EarthMoonOrbitPlane>({0 * SIUnit<Length>(),
                                             4E8 * SIUnit<Length>(),
                                             0 * SIUnit<Length>()}));
    Length const semi_major_axis = (q1 - q2).Norm();
    period_ = 2 * π * Sqrt(Pow<3>(semi_major_axis) /
                               (body1_->gravitational_parameter() +
                                body2_->gravitational_parameter()));
    centre_of_mass_ =
        Barycentre<Vector<Length, EarthMoonOrbitPlane>, Mass>(
            {q1, q2}, {body1_->mass(), body2_->mass()});
    Velocity<EarthMoonOrbitPlane> const v1(
        {-2 * π * (q1 - centre_of_mass_).Norm() / period_,
         0 * SIUnit<Speed>(),
         0 * SIUnit<Speed>()});
    Velocity<EarthMoonOrbitPlane> const v2(
        {2 * π * (q2 - centre_of_mass_).Norm() / period_,
         0 * SIUnit<Speed>(),
         0 * SIUnit<Speed>()});
    trajectory1_->Append(Instant(0 * SIUnit<Time>()), {q1, v1});
    trajectory2_->Append(Instant(0 * SIUnit<Time>()), {q2, v2});
    system_ = std::make_unique<NBodySystem<EarthMoonOrbitPlane>>();
  }

  template<typename Scalar, typename Frame>
  std::string ToMathematicaString(Vector<Scalar, Frame> const& vector) {
    R3Element<Scalar> const& coordinates = vector.coordinates();
    std::string result = "{";
    result += quantities::DebugString(coordinates.x);
    result += ",";
    result += quantities::DebugString(coordinates.y);
    result += ",";
    result += quantities::DebugString(coordinates.z);
    result += "}";
    return result;
  }

  template<typename Scalar, typename Frame>
  std::string ToMathematicaString(
      std::vector<Vector<Scalar, Frame>> const& vectors) {
    static std::string const mathematica_line =
        "\n(*****************************************************)\n";
    std::string result = mathematica_line;
    result += "ToExpression[StringReplace[\"\n{";
    std::string separator = "";
    for (const auto& vector : vectors) {
      result += separator;
      result += ToMathematicaString(vector);
      separator = ",\n";
    }
    result +=
        "}\",\n{\" m\"->\"\",\"e\"->\"*^\", \"\\n\"->\"\", \" \"->\"\"}]];";
    result += mathematica_line;
    return result;
  }

  template<typename T1, typename T2>
  std::vector<T2> ValuesOf(std::map<T1, Point<T2>> const& m,
                           Point<T2> const& relative_to) {
    std::vector<T2> result;
    for (auto const it : m) {
      result.push_back(it.second - relative_to);
    }
    return result;
  }

  std::unique_ptr<Body> body1_;
  std::unique_ptr<Body> body2_;
  std::unique_ptr<Body> body3_;
  std::unique_ptr<Trajectory<EarthMoonOrbitPlane>> trajectory1_;
  std::unique_ptr<Trajectory<EarthMoonOrbitPlane>> trajectory2_;
  std::unique_ptr<Trajectory<EarthMoonOrbitPlane>> trajectory3_;
  Position<EarthMoonOrbitPlane> centre_of_mass_;
  SPRKIntegrator<Length, Speed> integrator_;
  Time period_;
  std::unique_ptr<NBodySystem<EarthMoonOrbitPlane>> system_;
};

using NBodySystemDeathTest = NBodySystemTest;

TEST_F(NBodySystemDeathTest, IntegrateError) {
  EXPECT_DEATH({
    system_->Integrate(integrator_,
                       trajectory1_->last_time() + period_,
                       period_ / 100,
                       1,      // sampling_period
                       false,  // tmax_is_exact
                       {trajectory1_.get(),
                        trajectory2_.get(),
                        trajectory1_.get()});
  }, DeathMessage("Multiple trajectories"));
  EXPECT_DEATH({
    std::unique_ptr<Trajectory<EarthMoonOrbitPlane>> trajectory(
        new Trajectory<EarthMoonOrbitPlane>(*body2_));
    trajectory->Append(Instant(1 * SIUnit<Time>()),
                       {Position<EarthMoonOrbitPlane>(),
                        Velocity<EarthMoonOrbitPlane>()});
    system_->Integrate(integrator_,
                       trajectory1_->last_time() + period_,
                       period_ / 100,
                       1,      // sampling_period
                       false,  // tmax_is_exact
                       {trajectory1_.get(), trajectory.get()});
  }, DeathMessage("Inconsistent last time"));
}

// The canonical Earth-Moon system, tuned to produce circular orbits.
TEST_F(NBodySystemTest, EarthMoon) {
  std::vector<Vector<Length, EarthMoonOrbitPlane>> positions;
  system_->Integrate(integrator_,
                     trajectory1_->last_time() + period_,
                     period_ / 100,
                     1,      // sampling_period
                     false,  // tmax_is_exact
                     {trajectory1_.get(), trajectory2_.get()});

  positions = ValuesOf(trajectory1_->Positions(), centre_of_mass_);
  EXPECT_THAT(positions.size(), Eq(101));
  LOG(INFO) << ToMathematicaString(positions);
  EXPECT_THAT(Abs(positions[25].coordinates().y), Lt(3E-2 * SIUnit<Length>()));
  EXPECT_THAT(Abs(positions[50].coordinates().x), Lt(3E-2 * SIUnit<Length>()));
  EXPECT_THAT(Abs(positions[75].coordinates().y), Lt(3E-2 * SIUnit<Length>()));
  EXPECT_THAT(Abs(positions[100].coordinates().x), Lt(3E-2 * SIUnit<Length>()));

  positions = ValuesOf(trajectory2_->Positions(), centre_of_mass_);
  LOG(INFO) << ToMathematicaString(positions);
  EXPECT_THAT(positions.size(), Eq(101));
  EXPECT_THAT(Abs(positions[25].coordinates().y), Lt(2 * SIUnit<Length>()));
  EXPECT_THAT(Abs(positions[50].coordinates().x), Lt(2 * SIUnit<Length>()));
  EXPECT_THAT(Abs(positions[75].coordinates().y), Lt(2 * SIUnit<Length>()));
  EXPECT_THAT(Abs(positions[100].coordinates().x), Lt(2 * SIUnit<Length>()));
}

// Same as above, but the trajectories are passed in the reverse order.
TEST_F(NBodySystemTest, MoonEarth) {
  std::vector<Vector<Length, EarthMoonOrbitPlane>> positions;
  system_->Integrate(integrator_,
                     trajectory1_->last_time() + period_,
                     period_ / 100,
                     1,      // sampling_period
                     false,  // tmax_is_exact
                     {trajectory2_.get(), trajectory1_.get()});

  positions = ValuesOf(trajectory1_->Positions(), centre_of_mass_);
  EXPECT_THAT(positions.size(), Eq(101));
  LOG(INFO) << ToMathematicaString(positions);
  EXPECT_THAT(Abs(positions[25].coordinates().y), Lt(3E-2 * SIUnit<Length>()));
  EXPECT_THAT(Abs(positions[50].coordinates().x), Lt(3E-2 * SIUnit<Length>()));
  EXPECT_THAT(Abs(positions[75].coordinates().y), Lt(3E-2 * SIUnit<Length>()));
  EXPECT_THAT(Abs(positions[100].coordinates().x), Lt(3E-2 * SIUnit<Length>()));

  positions = ValuesOf(trajectory2_->Positions(), centre_of_mass_);
  LOG(INFO) << ToMathematicaString(positions);
  EXPECT_THAT(positions.size(), Eq(101));
  EXPECT_THAT(Abs(positions[25].coordinates().y), Lt(2 * SIUnit<Length>()));
  EXPECT_THAT(Abs(positions[50].coordinates().x), Lt(2 * SIUnit<Length>()));
  EXPECT_THAT(Abs(positions[75].coordinates().y), Lt(2 * SIUnit<Length>()));
  EXPECT_THAT(Abs(positions[100].coordinates().x), Lt(2 * SIUnit<Length>()));
}

// The Moon alone.  It moves in straight line.
TEST_F(NBodySystemTest, Moon) {
  Position<EarthMoonOrbitPlane> const reference_position;
  system_->Integrate(integrator_,
                     trajectory1_->last_time() + period_,
                     period_ / 100,
                     1,      // sampling_period
                     false,  // tmax_is_exact
                     {trajectory2_.get()});

  Length const q2 =
      (trajectory2_->last_position() - reference_position).coordinates().y;
  Speed const v2 = trajectory2_->last_velocity().coordinates().x;
  std::vector<Vector<Length, EarthMoonOrbitPlane>> const positions =
      ValuesOf(trajectory2_->Positions(), reference_position);
  LOG(INFO) << ToMathematicaString(positions);
  EXPECT_THAT(positions.size(), Eq(101));
  EXPECT_THAT(positions[25].coordinates().x, Eq(0.25 * period_ * v2));
  EXPECT_THAT(positions[25].coordinates().y, Eq(q2));
  EXPECT_THAT(positions[50].coordinates().x, Eq(0.50 * period_ * v2));
  EXPECT_THAT(positions[50].coordinates().y, Eq(q2));
  EXPECT_THAT(positions[75].coordinates().x,
              AlmostEquals(0.75 * period_ * v2, 1));
  EXPECT_THAT(positions[75].coordinates().y, Eq(q2));
  EXPECT_THAT(positions[100].coordinates().x, Eq(1.00 * period_ * v2));
  EXPECT_THAT(positions[100].coordinates().y, Eq(q2));
}

// The Earth and a massless probe 1 billion meters away, with the same velocity,
// and an acceleration which exactly compensates gravitational attraction.  Both
// bodies move in straight lines.
TEST_F(NBodySystemTest, EarthProbe) {
  Position<EarthMoonOrbitPlane> const reference_position;
  Length const distance = 1E9 * SIUnit<Length>();
  trajectory3_->Append(trajectory1_->last_time(),
                       {trajectory1_->last_position() +
                            Vector<Length, EarthMoonOrbitPlane>(
                                {0 * SIUnit<Length>(),
                                 distance,
                                 0 * SIUnit<Length>()}),
                        trajectory1_->last_velocity()});
  trajectory3_->set_intrinsic_acceleration(
      [this, distance](Instant const& t) {
    return Vector<Acceleration, EarthMoonOrbitPlane>(
        {0 * SIUnit<Acceleration>(),
         body1_->gravitational_parameter() / (distance * distance),
         0 * SIUnit<Acceleration>()});});

  system_->Integrate(integrator_,
                     trajectory1_->last_time() + period_,
                     period_ / 100,
                     1,      // sampling_period
                     false,  // tmax_is_exact
                     {trajectory1_.get(), trajectory3_.get()});

  Length const q1 =
      (trajectory1_->last_position() - reference_position).coordinates().y;
  Speed const v1 = trajectory1_->last_velocity().coordinates().x;
  std::vector<Vector<Length, EarthMoonOrbitPlane>> const positions1 =
      ValuesOf(trajectory1_->Positions(), reference_position);
  LOG(INFO) << ToMathematicaString(positions1);
  EXPECT_THAT(positions1.size(), Eq(101));
  EXPECT_THAT(positions1[25].coordinates().x,
              AlmostEquals(0.25 * period_ * v1, 1));
  EXPECT_THAT(positions1[25].coordinates().y, Eq(q1));
  EXPECT_THAT(positions1[50].coordinates().x,
              AlmostEquals(0.50 * period_ * v1, 1));
  EXPECT_THAT(positions1[50].coordinates().y, Eq(q1));
  EXPECT_THAT(positions1[75].coordinates().x,
              AlmostEquals(0.75 * period_ * v1, 1));
  EXPECT_THAT(positions1[75].coordinates().y, Eq(q1));
  EXPECT_THAT(positions1[100].coordinates().x,
              AlmostEquals(1.00 * period_ * v1, 1));
  EXPECT_THAT(positions1[100].coordinates().y, Eq(q1));

  Length const q3 =
      (trajectory3_->last_position() - reference_position).coordinates().y;
  Speed const v3 = trajectory3_->last_velocity().coordinates().x;
  std::vector<Vector<Length, EarthMoonOrbitPlane>> const positions3 =
      ValuesOf(trajectory3_->Positions(), reference_position);
  LOG(INFO) << ToMathematicaString(positions3);
  EXPECT_THAT(positions3.size(), Eq(101));
  EXPECT_THAT(positions3[25].coordinates().x,
              AlmostEquals(0.25 * period_ * v3, 1));
  EXPECT_THAT(positions3[25].coordinates().y, AlmostEquals(q3, 2));
  EXPECT_THAT(positions3[50].coordinates().x,
              AlmostEquals(0.50 * period_ * v3, 1));
  EXPECT_THAT(positions3[50].coordinates().y, AlmostEquals(q3, 2));
  EXPECT_THAT(positions3[75].coordinates().x,
              AlmostEquals(0.75 * period_ * v3, 1));
  EXPECT_THAT(positions3[75].coordinates().y, AlmostEquals(q3, 1));
  EXPECT_THAT(positions3[100].coordinates().x,
              AlmostEquals(1.00 * period_ * v3, 1));
  EXPECT_THAT(positions3[100].coordinates().y, Eq(q3));
}

TEST_F(NBodySystemTest, Sputnik1ToSputnik2) {
  std::unique_ptr<SolarSystem> const evolved_system =
      SolarSystem::AtСпутник1Launch(
          SolarSystem::Accuracy::kMinorAndMajorBodies);
  std::unique_ptr<SolarSystem> const at_спутник_2_launch =
      SolarSystem::AtСпутник2Launch(
          SolarSystem::Accuracy::kMinorAndMajorBodies);
  NBodySystem<ICRFJ2000Ecliptic> system;
  system.Integrate(
      integrator_,
      at_спутник_2_launch->trajectories().front()->last_time(),  // tmax
      45 * Minute,  // Δt
      0,  // sampling_period
      true,  // tmax_is_exact
      evolved_system->trajectories());  // trajectories
  for (std::size_t i = 0; i < evolved_system->trajectories().size(); ++i) {
    double const position_error = RelativeError(
        at_спутник_2_launch->trajectories()[i]->last_position() -
            kSolarSystemBarycentre,
        evolved_system->trajectories()[i]->last_position() -
            kSolarSystemBarycentre);
    double const velocity_error = RelativeError(
        at_спутник_2_launch->trajectories()[i]->last_velocity(),
        evolved_system->trajectories()[i]->last_velocity());
    // Upper bounds, tight to the nearest order of magnitude.
    double expected_velocity_error = 0;
    double expected_position_error = 0;
    Angle expected_angle_error = 0 * Degree;
    switch (i) {
      case SolarSystem::kSun:
        expected_position_error = 1E-8;
        expected_velocity_error = 1E-7;
        break;
      case SolarSystem::kJupiter:
        expected_position_error = 1E-8;
        expected_velocity_error = 1E-7;
        break;
      case SolarSystem::kSaturn:
        expected_position_error = 1E-8;
        expected_velocity_error = 1E-7;
        break;
      case SolarSystem::kNeptune:
        expected_position_error = 1E-8;
        expected_velocity_error = 1E-6;
        break;
      case SolarSystem::kUranus:
        expected_position_error = 1E-8;
        expected_velocity_error = 1E-5;
        break;
      case SolarSystem::kEarth:
        expected_position_error = 1E-8;
        expected_velocity_error = 1E-7;
        break;
      case SolarSystem::kVenus:
        expected_position_error = 1E-7;
        expected_velocity_error = 1E-7;
        break;
      case SolarSystem::kMars:
        expected_position_error = 1E-9;
        expected_velocity_error = 1E-8;
        break;
      case SolarSystem::kMercury:
        // NOTE(egg): We expect the error here to be higher than for the other
        // planets.
        expected_position_error = 1E-6;
        expected_velocity_error = 1E-6;
        break;
      case SolarSystem::kGanymede:
        expected_position_error = 1E-8;
        expected_velocity_error = 1E-4;
        expected_angle_error = 0.28 * Degree;
        break;
      case SolarSystem::kTitan:
        expected_position_error = 1E-5;
        expected_velocity_error = 1E-3;
        expected_angle_error = 0.089 * Degree;
        break;
      case SolarSystem::kCallisto:
        expected_position_error = 1E-8;
        expected_velocity_error = 1E-6;
        break;
      case SolarSystem::kIo:
        expected_position_error = 1E-6;
        expected_velocity_error = 1E-3;
        expected_angle_error = 7.6 * Degree;
        break;
      case SolarSystem::kMoon:
        expected_position_error = 1E-7;
        expected_velocity_error = 1E-6;
        break;
      case SolarSystem::kEuropa:
        expected_position_error = 1E-7;
        expected_velocity_error = 1E-3;
        expected_angle_error = 1.4 * Degree;
        break;
      case SolarSystem::kTriton:
        expected_position_error = 1E-7;
        expected_velocity_error = 1E-3;
        break;
      case SolarSystem::kEris:
        // NOTE(egg): we may want Dysnomia.
        expected_position_error = 1E-5;
        expected_velocity_error = 1E-5;
        break;
      case SolarSystem::kPluto:
        expected_position_error = 1E-6;
        expected_velocity_error = 1E-3;
        break;
      case SolarSystem::kTitania:
        expected_position_error = 1E-6;
        expected_velocity_error = 1E-3;
        break;
      case SolarSystem::kOberon:
        expected_position_error = 1E-8;
        expected_velocity_error = 1E-4;
        break;
      case SolarSystem::kRhea:
        expected_position_error = 1E-5;
        expected_velocity_error = 1E-1;
        expected_angle_error = 1.5 * Degree;
        break;
      case SolarSystem::kIapetus:
        expected_position_error = 1E-7;
        expected_velocity_error = 1E-5;
        break;
      case SolarSystem::kCharon:
        expected_position_error = 1E-6;
        expected_velocity_error = 1E-3;
        break;
      case SolarSystem::kAriel:
        expected_position_error = 1E-6;
        expected_velocity_error = 1E-2;
        expected_angle_error = 0.74 * Degree;
        break;
      case SolarSystem::kUmbriel:
        expected_position_error = 1E-6;
        expected_velocity_error = 1E-2;
        expected_angle_error = 0.24 * Degree;
        break;
      case SolarSystem::kDione:
        expected_position_error = 1E-4;
        expected_velocity_error = 1E-0;
        expected_angle_error = 4.8 * Degree;
        break;
      case SolarSystem::kTethys:
        expected_position_error = 1E-4;
        expected_velocity_error = 1E-0;
        expected_angle_error = 11.4 * Degree;
        break;
      default:
        LOG(FATAL) << "No such body";
    }
    EXPECT_THAT(position_error, Lt(expected_position_error))
        << SolarSystem::name(i);
    EXPECT_THAT(position_error, Gt(expected_position_error / 10.0))
        << SolarSystem::name(i);
    EXPECT_THAT(velocity_error, Lt(expected_velocity_error))
        << SolarSystem::name(i);
    EXPECT_THAT(velocity_error, Gt(expected_velocity_error / 10.0))
        << SolarSystem::name(i);
    if (i != SolarSystem::kSun) {
      // Look at the error in the position relative to the parent.
      Vector<Length, ICRFJ2000Ecliptic> expected =
          at_спутник_2_launch->trajectories()[i]->last_position() -
              at_спутник_2_launch->trajectories()[SolarSystem::parent(i)]->
                  last_position();
      Vector<Length, ICRFJ2000Ecliptic> actual =
          evolved_system->trajectories()[i]->last_position() -
              evolved_system->trajectories()[SolarSystem::parent(i)]->
                  last_position();
      double const vector_error =  RelativeError(expected, actual);
      if (i == SolarSystem::kIo ||
          i == SolarSystem::kEuropa ||
          i == SolarSystem::kGanymede ||
          i == SolarSystem::kTitan ||
          i == SolarSystem::kRhea ||
          i == SolarSystem::kAriel ||
          i == SolarSystem::kUmbriel ||
          i == SolarSystem::kDione ||
          i == SolarSystem::kTethys) {
        double const length_error = RelativeError(expected.Norm(),
                                                  actual.Norm());
        Area const product_of_norms = expected.Norm() * actual.Norm();
        Angle const angle = ArcTan(
            Wedge(expected, actual).Norm() / product_of_norms,
            InnerProduct(expected, actual) / product_of_norms);
        // We are missing some sort of precession here, probably due to
        // quadrupole moment.
        LOG(ERROR) << SolarSystem::name(SolarSystem::parent(i));
        LOG(ERROR) << angle * Pow<4>(expected.Norm());
        EXPECT_THAT(angle / Degree,
                    Gt(expected_angle_error / Degree * 0.9))
            << SolarSystem::name(i);
        EXPECT_THAT(angle / Degree,
                    Lt(expected_angle_error / Degree * 1.1))
            << SolarSystem::name(i);
        EXPECT_THAT(length_error, Lt(3E-3)) << SolarSystem::name(i);
        EXPECT_THAT(length_error, Gt(1E-6)) << SolarSystem::name(i);
      } else if (SolarSystem::parent(i) != SolarSystem::kSun &&
                 SolarSystem::parent(i) != SolarSystem::kEarth) {
        // Less extreme precession, probably the same cause.
        EXPECT_THAT(vector_error, Lt(1E-3)) << SolarSystem::name(i);
        EXPECT_THAT(vector_error, Gt(1E-5)) << SolarSystem::name(i);
      } else if (SolarSystem::parent(i) == SolarSystem::kPluto) {
        // We are missing Hydra and Nyx.
        EXPECT_THAT(vector_error, Lt(1E-4)) << SolarSystem::name(i);
        EXPECT_THAT(vector_error, Gt(1E-5)) << SolarSystem::name(i);
      } else if (i == SolarSystem::kMoon) {
        // What is this?
        EXPECT_THAT(vector_error, Lt(1E-5)) << SolarSystem::name(i);
        EXPECT_THAT(vector_error, Gt(1E-6)) << SolarSystem::name(i);
      } else if (i == SolarSystem::kEris) {
        // We are missing Dysnomia.
        EXPECT_THAT(vector_error, Lt(1E-5)) << SolarSystem::name(i);
        EXPECT_THAT(vector_error, Gt(1E-6)) << SolarSystem::name(i);
      } else if (i == SolarSystem::kPluto) {
        // We are missing Hydra and Nyx.
        EXPECT_THAT(vector_error, Lt(1E-6)) << SolarSystem::name(i);
        EXPECT_THAT(vector_error, Gt(1E-7)) << SolarSystem::name(i);
      } else if (i == SolarSystem::kMercury) {
        // We are missing GR...
        EXPECT_THAT(vector_error, Lt(1E-6)) << SolarSystem::name(i);
        EXPECT_THAT(vector_error, Gt(1E-7)) << SolarSystem::name(i);
      } else {
        EXPECT_THAT(vector_error, Lt(1E-7)) << SolarSystem::name(i);
      }
    }
  }
}

}  // namespace physics
}  // namespace principia
