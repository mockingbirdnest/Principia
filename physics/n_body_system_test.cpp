#include "physics/n_body_system.hpp"

#include <memory>
#include <string>
#include <vector>

#include "geometry/grassmann.hpp"
#include "geometry/point.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "physics/body.hpp"
#include "quantities/constants.hpp"
#include "quantities/numbers.hpp"

using principia::constants::GravitationalConstant;
using principia::geometry::Barycentre;
using principia::geometry::Point;
using principia::geometry::Vector;
using principia::quantities::Pow;
using principia::quantities::SIUnit;
using testing::Eq;
using testing::Lt;

namespace principia {
namespace physics {

class NBodySystemTest : public testing::Test {
 protected:
  struct EarthMoonBarycentricFrame;

  void SetUp() override {
    integrator_.Initialize(integrator_.Order5Optimal());

    // The Earth-Moon system, roughly, with a circular orbit with velocities
    // in the center-of-mass frame.
    body1_ = new Body<EarthMoonBarycentricFrame>(6E24 * SIUnit<Mass>());
    body2_ = new Body<EarthMoonBarycentricFrame>(7E22 * SIUnit<Mass>());
    Point<Vector<Length, EarthMoonBarycentricFrame>> const
        q1(Vector<Length, EarthMoonBarycentricFrame>({0 * SIUnit<Length>(),
                                          0 * SIUnit<Length>(),
                                          0 * SIUnit<Length>()}));
    Point<Vector<Length, EarthMoonBarycentricFrame>> const
        q2(Vector<Length, EarthMoonBarycentricFrame>({0 * SIUnit<Length>(),
                                          4E8 * SIUnit<Length>(),
                                          0 * SIUnit<Length>()}));
    Point<Vector<Length, EarthMoonBarycentricFrame>> const centre_of_mass =
        Barycentre(q1, body1_->mass(), q2, body2_->mass());
    Length const semi_major_axis = (q1 - q2).Norm();
    period_ = 2 * π * Sqrt(Pow<3>(semi_major_axis) /
                               (body1_->gravitational_parameter() +
                                body2_->gravitational_parameter()));
    Point<Vector<Speed, EarthMoonBarycentricFrame>> const
        v1(Vector<Speed, EarthMoonBarycentricFrame>({
            -2 * π * (q1 - centre_of_mass).Norm() / period_,
            0 * SIUnit<Speed>(),
            0 * SIUnit<Speed>()}));
    Point<Vector<Speed, EarthMoonBarycentricFrame>> const
        v2(Vector<Speed, EarthMoonBarycentricFrame>({
            2 * π * (q2 - centre_of_mass).Norm() / period_,
            0 * SIUnit<Speed>(),
            0 * SIUnit<Speed>()}));
    Point<Vector<Speed, EarthMoonBarycentricFrame>> const overall_velocity =
        Barycentre(v1, body1_->mass(), v2, body2_->mass());
    body1_->AppendToTrajectory(q1 - centre_of_mass,
                               v1 - overall_velocity,
                               0 * SIUnit<Time>());
    body2_->AppendToTrajectory(q2 - centre_of_mass,
                               v2 - overall_velocity,
                               0 * SIUnit<Time>());
    NBodySystem<EarthMoonBarycentricFrame>::Bodies* massive_bodies = 
        new NBodySystem<EarthMoonBarycentricFrame>::Bodies;
    massive_bodies->emplace_back(body1_);
    massive_bodies->emplace_back(body2_);
    system_.reset(new NBodySystem<EarthMoonBarycentricFrame>(massive_bodies,
                                                             nullptr));
  }

  template<typename Scalar, typename Frame>
  std::string ToMathematicaString(Vector<Scalar, Frame> const& vector) {
    R3Element<Scalar> const& coordinates = vector.coordinates();
    std::string result = "{";
    result += quantities::ToString(coordinates.x);
    result += ",";
    result += quantities::ToString(coordinates.y);
    result += ",";
    result += quantities::ToString(coordinates.z);
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

  Body<EarthMoonBarycentricFrame>* body1_;
  Body<EarthMoonBarycentricFrame>* body2_;
  SPRKIntegrator<Length, Speed> integrator_;
  Time period_;
  std::unique_ptr<NBodySystem<EarthMoonBarycentricFrame>> system_;
};

TEST_F(NBodySystemTest, EarthMoon) {
  std::vector<Vector<Length, EarthMoonBarycentricFrame>> positions;
  system_->Integrate(integrator_, period_, period_ / 100, 1);

  positions = body1_->positions();
  EXPECT_THAT(positions.size(), Eq(101));
  LOG(INFO) << ToMathematicaString(positions);
  EXPECT_THAT(Abs(positions[25].coordinates().y), Lt(3E-2 * SIUnit<Length>()));
  EXPECT_THAT(Abs(positions[50].coordinates().x), Lt(3E-2 * SIUnit<Length>()));
  EXPECT_THAT(Abs(positions[75].coordinates().y), Lt(3E-2 * SIUnit<Length>()));
  EXPECT_THAT(Abs(positions[100].coordinates().x), Lt(3E-2 * SIUnit<Length>()));

  positions = body2_->positions();
  LOG(INFO) << ToMathematicaString(positions);
  EXPECT_THAT(positions.size(), Eq(101));
  EXPECT_THAT(Abs(positions[25].coordinates().y), Lt(2 * SIUnit<Length>()));
  EXPECT_THAT(Abs(positions[50].coordinates().x), Lt(2 * SIUnit<Length>()));
  EXPECT_THAT(Abs(positions[75].coordinates().y), Lt(2 * SIUnit<Length>()));
  EXPECT_THAT(Abs(positions[100].coordinates().x), Lt(2 * SIUnit<Length>()));
}

}  // namespace physics
}  // namespace principia
