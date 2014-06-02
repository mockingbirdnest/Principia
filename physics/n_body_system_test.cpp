#include "physics/n_body_system.hpp"

#include <memory>
#include <string>

#include "geometry/grassmann.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "physics/body.hpp"
#include "physics/frame.hpp"

using principia::geometry::Vector;

namespace principia {
namespace physics {

class NBodySystemTest : public testing::Test {
 protected:
  void SetUp() override {
    integrator_.Initialize(integrator_.Order5Optimal());

    //body1_ = new Body<InertialFrame>(1 * GravitationalParameter::SIUnit());
    //body2_ = new Body<InertialFrame>(2 * GravitationalParameter::SIUnit());
    //Vector<Length, InertialFrame> q1({1 * Length::SIUnit(),
    //                                  2 * Length::SIUnit(),
    //                                  3 * Length::SIUnit()});
    //Vector<Length, InertialFrame> q2({4 * Length::SIUnit(),
    //                                  5 * Length::SIUnit(),
    //                                  6 * Length::SIUnit()});
    //Vector<Momentum, InertialFrame> p1({3 * Momentum::SIUnit(),
    //                                    2 * Momentum::SIUnit(),
    //                                    1 * Momentum::SIUnit()});
    //Vector<Momentum, InertialFrame> p2({6 * Momentum::SIUnit(),
    //                                    5 * Momentum::SIUnit(),
    //                                    4 * Momentum::SIUnit()});
    body1_ = new Body<InertialFrame>(1 * GravitationalParameter::SIUnit());
    body2_ = new Body<InertialFrame>(2 * GravitationalParameter::SIUnit());
    Vector<Length, InertialFrame> q1({1 * Length::SIUnit(),
                                      0 * Length::SIUnit(),
                                      0 * Length::SIUnit()});
    Vector<Length, InertialFrame> q2({0 * Length::SIUnit(),
                                      2 * Length::SIUnit(),
                                      0 * Length::SIUnit()});
    Vector<Momentum, InertialFrame> p1({0 * Momentum::SIUnit(),
                                        2 * Momentum::SIUnit(),
                                        0 * Momentum::SIUnit()});
    Vector<Momentum, InertialFrame> p2({1 * Momentum::SIUnit(),
                                        0 * Momentum::SIUnit(),
                                        0 * Momentum::SIUnit()});
    body1_->AppendToTrajectory({q1}, {p1}, {0 * Time::SIUnit()});
    body2_->AppendToTrajectory({q2}, {p2}, {0 * Time::SIUnit()});
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
  std::unique_ptr<NBodySystem> system_;
};

TEST_F(NBodySystemTest, T) {
  std::vector<Vector<Length, InertialFrame>> positions;
  std::vector<Vector<Momentum, InertialFrame>> momenta;
  std::vector<Time> times;
  system_->Integrate(integrator_, 10 * Time::SIUnit(), 0.1 * Time::SIUnit(), 1);
  body1_->GetTrajectory(&positions, &momenta, &times);
  LOG(ERROR) << ToMathematicaString(positions);
  body2_->GetTrajectory(&positions, &momenta, &times);
  LOG(ERROR) << ToMathematicaString(positions);
}

}  // namespace physics
}  // namespace principia