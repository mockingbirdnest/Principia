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

    // The Earth-Moon system, roughly.
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
    Point<Vector<Speed, InertialFrame>> const
        v1(Vector<Speed, InertialFrame>({0 * Speed::SIUnit(),
                                         0 * Speed::SIUnit(),
                                         0 * Speed::SIUnit()}));
    // 1E3 m/s puts us at the apogee at the beginning, because the speed for a
    // circular orbit with this separation would be 1006.36 m/s.
    Point<Vector<Speed, InertialFrame>> const
        v2(Vector<Speed, InertialFrame>({1E3 * Speed::SIUnit(),
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
    Length const r2 = (q2 - centre_of_mass).Norm();
    Length const semi_major_axis =
        1 / (2 / r2 -
                 (v2 - overall_velocity).Norm().Pow<2>() /
                      (body1_->gravitational_parameter() +
                       body2_->gravitational_parameter()));
    Length const r1 = (q1 - centre_of_mass).Norm();
    Length const semi_major_axis1 =
        1 / (2 / r1 -
                 (v1 - overall_velocity).Norm().Pow<2>() /
                      (body1_->gravitational_parameter() +
                       body2_->gravitational_parameter()));
    LOG(ERROR)<<semi_major_axis<<" "<<semi_major_axis1;
    period_ = 2 * π * Sqrt(
        (semi_major_axis1 + semi_major_axis).Pow<3>() /
            (GravitationalConstant * (body1_->mass() + body2_->mass())));
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
  LOG(ERROR)<<period_;
  system_->Integrate(integrator_, 1.1*period_, period_ / 1000, 1);
  body1_->GetTrajectory(&positions, &momenta, &times);
  LOG(ERROR) << ToMathematicaString(positions);
  body2_->GetTrajectory(&positions, &momenta, &times);
  LOG(ERROR) << ToMathematicaString(positions);
}

}  // namespace physics
}  // namespace principia