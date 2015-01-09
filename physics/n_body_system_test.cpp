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
#include "physics/massive_body.hpp"
#include "physics/massless_body.hpp"
#include "physics/trajectory.hpp"
#include "quantities/constants.hpp"
#include "quantities/numbers.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/almost_equals.hpp"
#include "testing_utilities/numerics.hpp"
#include "testing_utilities/solar_system.hpp"

using principia::constants::GravitationalConstant;
using principia::geometry::Instant;
using principia::geometry::Point;
using principia::geometry::Vector;
using principia::quantities::Angle;
using principia::quantities::ArcTan;
using principia::quantities::Area;
using principia::quantities::Pow;
using principia::quantities::SIUnit;
using principia::testing_utilities::AlmostEquals;
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
  enum class Tag {
    kEarthMoonOrbitPlane,
  };

  using EarthMoonOrbitPlane = Frame<Tag, Tag::kEarthMoonOrbitPlane, true>;

  NBodySystemTest()
      : body1_(MassiveBody(6E24 * SIUnit<Mass>())),
        body2_(MassiveBody(7E22 * SIUnit<Mass>())) {
    integrator_.Initialize(integrator_.Order5Optimal());

    // The Earth-Moon system, roughly, with a circular orbit with velocities
    // in the centre-of-mass frame.
    trajectory1_ = std::make_unique<Trajectory<EarthMoonOrbitPlane>>(
                       check_not_null(&body1_));
    trajectory2_ = std::make_unique<Trajectory<EarthMoonOrbitPlane>>(
                       check_not_null(&body2_));
    trajectory3_ = std::make_unique<Trajectory<EarthMoonOrbitPlane>>(
                       check_not_null(&body3_));
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
                               (body1_.gravitational_parameter() +
                                body2_.gravitational_parameter()));
    centre_of_mass_ =
        geometry::Barycentre<Vector<Length, EarthMoonOrbitPlane>, Mass>(
            {q1, q2}, {body1_.mass(), body2_.mass()});
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
    for (auto const& vector : vectors) {
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

  MassiveBody body1_;
  MassiveBody body2_;
  MasslessBody body3_;  // A massless probe.
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
                       trajectory1_->last().time() + period_,
                       period_ / 100,
                       1,      // sampling_period
                       false,  // tmax_is_exact
                       {trajectory1_.get(),
                        trajectory2_.get(),
                        trajectory1_.get()});
  }, "Multiple trajectories");
  EXPECT_DEATH({
    std::unique_ptr<Trajectory<EarthMoonOrbitPlane>> trajectory =
        std::make_unique<Trajectory<EarthMoonOrbitPlane>>(check_not_null(&body2_));
    trajectory->Append(Instant(1 * SIUnit<Time>()),
                       {Position<EarthMoonOrbitPlane>(),
                        Velocity<EarthMoonOrbitPlane>()});
    system_->Integrate(integrator_,
                       trajectory1_->last().time() + period_,
                       period_ / 100,
                       1,      // sampling_period
                       false,  // tmax_is_exact
                       {trajectory1_.get(), trajectory.get()});
  }, "Inconsistent last time");
}

// The canonical Earth-Moon system, tuned to produce circular orbits.
TEST_F(NBodySystemTest, EarthMoon) {
  std::vector<Vector<Length, EarthMoonOrbitPlane>> positions;
  system_->Integrate(integrator_,
                     trajectory1_->last().time() + period_,
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
                     trajectory1_->last().time() + period_,
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
  Position<EarthMoonOrbitPlane> const reference_position =
      Position<EarthMoonOrbitPlane>();
  system_->Integrate(integrator_,
                     trajectory1_->last().time() + period_,
                     period_ / 100,
                     1,      // sampling_period
                     false,  // tmax_is_exact
                     {trajectory2_.get()});

  Length const q2 = (trajectory2_->last().degrees_of_freedom().position() -
                     reference_position).coordinates().y;
  Speed const v2 =
      trajectory2_->last().degrees_of_freedom().velocity().coordinates().x;
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
  Position<EarthMoonOrbitPlane> const reference_position =
      Position<EarthMoonOrbitPlane>();
  Length const distance = 1E9 * SIUnit<Length>();
  trajectory3_->Append(trajectory1_->last().time(),
                       {trajectory1_->last().degrees_of_freedom().position() +
                            Vector<Length, EarthMoonOrbitPlane>(
                                {0 * SIUnit<Length>(),
                                 distance,
                                 0 * SIUnit<Length>()}),
                        trajectory1_->last().degrees_of_freedom().velocity()});
  trajectory3_->set_intrinsic_acceleration(
      [this, distance](Instant const& t) {
    return Vector<Acceleration, EarthMoonOrbitPlane>(
        {0 * SIUnit<Acceleration>(),
         body1_.gravitational_parameter() / (distance * distance),
         0 * SIUnit<Acceleration>()});});

  system_->Integrate(integrator_,
                     trajectory1_->last().time() + period_,
                     period_ / 100,
                     1,      // sampling_period
                     false,  // tmax_is_exact
                     {trajectory1_.get(), trajectory3_.get()});

  Length const q1 = (trajectory1_->last().degrees_of_freedom().position() -
                     reference_position).coordinates().y;
  Speed const v1 =
      trajectory1_->last().degrees_of_freedom().velocity().coordinates().x;
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

  Length const q3 = (trajectory3_->last().degrees_of_freedom().position() -
                     reference_position).coordinates().y;
  Speed const v3 =
      trajectory3_->last().degrees_of_freedom().velocity().coordinates().x;
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
          SolarSystem::Accuracy::kAllBodiesAndOblateness);
  std::unique_ptr<SolarSystem> const at_спутник_2_launch =
      SolarSystem::AtСпутник2Launch(
          SolarSystem::Accuracy::kAllBodiesAndOblateness);
  NBodySystem<ICRFJ2000Ecliptic> system;
  system.Integrate(
      integrator_,
      at_спутник_2_launch->trajectories().front()->last().time(),  // tmax
      45 * Minute,  // Δt
      0,  // sampling_period
      true,  // tmax_is_exact
      evolved_system->trajectories());  // trajectories

  // Upper bounds, tight to the nearest order of magnitude.
  static std::map<SolarSystem::Index, Angle> const expected_angle_error = {{}};
  static std::map<SolarSystem::Index,
                  double> const expected_parent_distance_error = {{}};
  static std::map<SolarSystem::Index,
                  double> const expected_parent_offset_error = {
      {SolarSystem::kAriel, 1E-3},
      {SolarSystem::kDione, 1E-3},
      {SolarSystem::kIo, 1E-3},
      {SolarSystem::kOberon, 1E-3},
      {SolarSystem::kTethys, 1E-3},
      {SolarSystem::kTitania, 1E-3},
      {SolarSystem::kTriton, 1E-4},
      {SolarSystem::kCharon, 1E-4},
      {SolarSystem::kEuropa, 1E-4},
      {SolarSystem::kRhea, 1E-4},
      {SolarSystem::kTitan, 1E-4},
      {SolarSystem::kUmbriel, 1E-4},
      {SolarSystem::kEris, 1E-5},  // NOTE(egg): we may want Dysnomia.
      {SolarSystem::kGanymede, 1E-5},
      {SolarSystem::kIapetus, 1E-5},
      {SolarSystem::kMoon, 1E-5},  // What is this?
      {SolarSystem::kCallisto, 1E-6},
      {SolarSystem::kMercury, 1E-6},  // NOTE(egg): General relativity.
      {SolarSystem::kPluto, 1E-6},  // NOTE(egg): We are missing Hydra and Nyx.
      {SolarSystem::kVenus, 1E-7},
      {SolarSystem::kEarth, 1E-8},
      {SolarSystem::kJupiter, 1E-8},
      {SolarSystem::kNeptune, 1E-8},
      {SolarSystem::kSaturn, 1E-8},
      {SolarSystem::kUranus, 1E-8},
      {SolarSystem::kMars, 1E-9}};
  static std::map<SolarSystem::Index, double> const expected_position_error = {
      {SolarSystem::kEris, 1E-5},  // NOTE(egg): we may want Dysnomia.
      {SolarSystem::kCharon, 1E-6},
      {SolarSystem::kMercury, 1E-6},  // NOTE(egg): General relativity.
      {SolarSystem::kPluto, 1E-6},
      {SolarSystem::kTethys, 1E-6},
      {SolarSystem::kAriel, 1E-7},
      {SolarSystem::kDione, 1E-7},
      {SolarSystem::kEuropa, 1E-7},
      {SolarSystem::kIo, 1E-7},
      {SolarSystem::kMoon, 1E-7},
      {SolarSystem::kOberon, 1E-7},
      {SolarSystem::kRhea, 1E-7},
      {SolarSystem::kTitan, 1E-7},
      {SolarSystem::kTitania, 1E-7},
      {SolarSystem::kVenus, 1E-7},
      {SolarSystem::kCallisto, 1E-8},
      {SolarSystem::kEarth, 1E-8},
      {SolarSystem::kGanymede, 1E-8},
      {SolarSystem::kIapetus, 1E-8},
      {SolarSystem::kJupiter, 1E-8},
      {SolarSystem::kNeptune, 1E-8},
      {SolarSystem::kSaturn, 1E-8},
      {SolarSystem::kSun, 1E-8},
      {SolarSystem::kTriton, 1E-8},
      {SolarSystem::kUmbriel, 1E-8},
      {SolarSystem::kUranus, 1E-8},
      {SolarSystem::kMars, 1E-9}};
  static std::map<SolarSystem::Index, double> const expected_velocity_error = {
      {SolarSystem::kAriel, 1E-3},
      {SolarSystem::kCharon, 1E-3},
      {SolarSystem::kDione, 1E-3},
      {SolarSystem::kIo, 1E-3},
      {SolarSystem::kPluto, 1E-3},
      {SolarSystem::kTethys, 1E-3},
      {SolarSystem::kEuropa, 1E-4},
      {SolarSystem::kOberon, 1E-4},
      {SolarSystem::kRhea, 1E-4},
      {SolarSystem::kTitania, 1E-4},
      {SolarSystem::kTriton, 1E-4},
      {SolarSystem::kUmbriel, 1E-4},
      {SolarSystem::kEris, 1E-5},  // NOTE(egg): we may want Dysnomia.
      {SolarSystem::kGanymede, 1E-5},
      {SolarSystem::kTitan, 1E-5},
      {SolarSystem::kUranus, 1E-5},
      {SolarSystem::kCallisto, 1E-6},
      {SolarSystem::kIapetus, 1E-6},
      {SolarSystem::kMercury, 1E-6},  // NOTE(egg): General relativity.
      {SolarSystem::kMoon, 1E-6},
      {SolarSystem::kSaturn, 1E-6},
      {SolarSystem::kEarth, 1E-7},
      {SolarSystem::kJupiter, 1E-7},
      {SolarSystem::kNeptune, 1E-7},
      {SolarSystem::kSun, 1E-7},
      {SolarSystem::kVenus, 1E-7},
      {SolarSystem::kMars, 1E-8}};

  for (std::size_t i = 0; i < evolved_system->trajectories().size(); ++i) {
    SolarSystem::Index const index = static_cast<SolarSystem::Index>(i);
    double const position_error = RelativeError(
        at_спутник_2_launch->trajectories()[i]->
            last().degrees_of_freedom().position() - kSolarSystemBarycentre,
        evolved_system->trajectories()[i]->
            last().degrees_of_freedom().position() - kSolarSystemBarycentre);
    double const velocity_error = RelativeError(
        at_спутник_2_launch->trajectories()[i]->
            last().degrees_of_freedom().velocity(),
        evolved_system->trajectories()[i]->
            last().degrees_of_freedom().velocity());
    EXPECT_THAT(position_error, Lt(expected_position_error.at(index)))
        << SolarSystem::name(i);
    EXPECT_THAT(position_error, Gt(expected_position_error.at(index) / 10.0))
        << SolarSystem::name(i);
    EXPECT_THAT(velocity_error, Lt(expected_velocity_error.at(index)))
        << SolarSystem::name(i);
    EXPECT_THAT(velocity_error, Gt(expected_velocity_error.at(index) / 10.0))
        << SolarSystem::name(i);
    if (i != SolarSystem::kSun) {
      // Look at the error in the position relative to the parent.
      Vector<Length, ICRFJ2000Ecliptic> expected =
          at_спутник_2_launch->trajectories()[i]->
              last().degrees_of_freedom().position() -
          at_спутник_2_launch->trajectories()[SolarSystem::parent(i)]->
              last().degrees_of_freedom().position();
      Vector<Length, ICRFJ2000Ecliptic> actual =
          evolved_system->trajectories()[i]->
              last().degrees_of_freedom().position() -
          evolved_system->trajectories()[SolarSystem::parent(i)]->
              last().degrees_of_freedom().position();
      if (expected_angle_error.find(index) != expected_angle_error.end()) {
        Area const product_of_norms = expected.Norm() * actual.Norm();
        Angle const angle = ArcTan(
            Wedge(expected, actual).Norm() / product_of_norms,
            InnerProduct(expected, actual) / product_of_norms);
        EXPECT_THAT(angle / Degree,
                    Gt(expected_angle_error.at(index) / Degree * 0.9))
            << SolarSystem::name(i);
        EXPECT_THAT(angle / Degree,
                    Lt(expected_angle_error.at(index) / Degree * 1.1))
            << SolarSystem::name(i);
      }
      if (expected_parent_distance_error.find(index) !=
          expected_parent_distance_error.end()) {
        double const parent_distance_error = RelativeError(expected.Norm(),
                                                  actual.Norm());
        EXPECT_THAT(parent_distance_error,
                    Lt(expected_parent_distance_error.at(index)))
            << SolarSystem::name(i);
        EXPECT_THAT(parent_distance_error,
                    Gt(expected_parent_distance_error.at(index) / 10.0))
            << SolarSystem::name(i);
      }
      if (expected_parent_offset_error.find(index) !=
          expected_parent_offset_error.end()) {
        double const parent_offset_error =  RelativeError(expected, actual);
        EXPECT_THAT(parent_offset_error,
                    Lt(expected_parent_offset_error.at(index)))
            << SolarSystem::name(i);
        EXPECT_THAT(parent_offset_error,
                    Gt(expected_parent_offset_error.at(index) / 10.0))
            << SolarSystem::name(i);
      }
    }
  }
}

}  // namespace physics
}  // namespace principia
