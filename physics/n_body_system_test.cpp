#include "physics/n_body_system.hpp"

#include <memory>
#include <string>

#include "geometry/grassmann.hpp"
#include "geometry/point.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "physics/body.hpp"
#include "physics/frame.hpp"
#include "quantities/constants.hpp"
#include "quantities/numbers.hpp"

using principia::constants::GravitationalConstant;
using principia::geometry::Barycentre;
using principia::geometry::Point;
using principia::geometry::Vector;

namespace principia {
namespace physics {

class NBodySystemTest : public testing::Test {
 protected:
  void SetUp() override {
    integrator_.Initialize(integrator_.Order5Optimal());

    // The Earth-Moon system, roughly, with a circular orbit with velocities
    // in the center-of-mass frame.
    body1_ = new Body<InertialFrame>(6E24 * Mass::SIUnit());
    body2_ = new Body<InertialFrame>(7E22 * Mass::SIUnit());
    Point<Vector<Length, InertialFrame>> const
        q1(Vector<Length, InertialFrame>({0 * Length::SIUnit(),
                                          0 * Length::SIUnit(),
                                          0 * Length::SIUnit()}));
    Point<Vector<Length, InertialFrame>> const
        q2(Vector<Length, InertialFrame>({0 * Length::SIUnit(),
                                          4E8 * Length::SIUnit(),
                                          0 * Length::SIUnit()}));
    Point<Vector<Length, InertialFrame>> const centre_of_mass =
        Barycentre(q1, body1_->mass(), q2, body2_->mass());
    Length const semi_major_axis = (q1 - q2).Norm();
    period_ = 2 * π * Sqrt(semi_major_axis.Pow<3>() /
                               (body1_->gravitational_parameter() +
                                body2_->gravitational_parameter()));
    Point<Vector<Speed, InertialFrame>> const
        v1(Vector<Speed, InertialFrame>({
            -2 * π * (q1 - centre_of_mass).Norm() / period_,
            0 * Speed::SIUnit(),
            0 * Speed::SIUnit()}));
    Point<Vector<Speed, InertialFrame>> const
        v2(Vector<Speed, InertialFrame>({
            2 * π * (q2 - centre_of_mass).Norm() / period_,
            0 * Speed::SIUnit(),
            0 * Speed::SIUnit()}));
    Point<Vector<Speed, InertialFrame>> const overall_velocity =
        Barycentre(v1, body1_->mass(), v2, body2_->mass());
    body1_->AppendToTrajectory({q1 - centre_of_mass},
                               {v1 - overall_velocity},
                               {0 * Time::SIUnit()});
    body2_->AppendToTrajectory({q2 - centre_of_mass},
                               {v2 - overall_velocity},
                               {0 * Time::SIUnit()});
    system_.reset(new NBodySystem(
        new std::vector<Body<InertialFrame>*>({body1_, body2_})));
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
        "(*****************************************************)";
    std::string result = mathematica_line + "\n";
    result += "ToExpression[StringReplace[\"\n{";
    std::string separator = "";
    for (const auto& vector : vectors) {
      result += separator;
      result += ToMathematicaString(vector);
      separator = ",\n";
    }
    result +=
        "}\",\n{\" m\"->\"\",\"e\"->\"*^\", \"\\n\"->\"\", \" \"->\"\"}]];\n";
    result += mathematica_line;
    return result;
  }

  Body<InertialFrame>* body1_;
  Body<InertialFrame>* body2_;
  SPRKIntegrator integrator_;
  Time period_;
  std::unique_ptr<NBodySystem> system_;
};

TEST_F(NBodySystemTest, T) {
  std::vector<Vector<Length, InertialFrame>> positions;
  std::vector<Vector<Speed, InertialFrame>> momenta;
  std::vector<Time> times;
  system_->Integrate(integrator_, period_, period_ / 100, 1);
  body1_->GetTrajectory(&positions, &momenta, &times);
  LOG(ERROR) << ToMathematicaString(positions);
  body2_->GetTrajectory(&positions, &momenta, &times);
  LOG(ERROR) << ToMathematicaString(positions);
}

}  // namespace physics
}  // namespace principia