#include "physics/n_body_system.hpp"

#include <memory>

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "physics/body.hpp"
#include "physics/frame.hpp"

namespace principia {
namespace physics {

class NBodySystemTest : public testing::Test {
 protected:
  void SetUp() override {
    integrator_.Initialize(integrator_.Order5Optimal());

    Body<InertialFrame>* body1(
        new Body<InertialFrame>(1E-11 * GravitationalParameter::SIUnit()));
    Body<InertialFrame>* body2(
        new Body<InertialFrame>(2E-11 * GravitationalParameter::SIUnit()));
    Vector<Length, InertialFrame> q1({1 * Length::SIUnit(),
                                      2 * Length::SIUnit(),
                                      3 * Length::SIUnit()});
    Vector<Length, InertialFrame> q2({4 * Length::SIUnit(),
                                      5 * Length::SIUnit(),
                                      6 * Length::SIUnit()});
    Vector<Momentum, InertialFrame> p1({3 * Momentum::SIUnit(),
                                        2 * Momentum::SIUnit(),
                                        1 * Momentum::SIUnit()});
    Vector<Momentum, InertialFrame> p2({6 * Momentum::SIUnit(),
                                        5 * Momentum::SIUnit(),
                                        4 * Momentum::SIUnit()});
    body1->AppendToTrajectory({q1}, {p1}, {0 * Time::SIUnit()});
    body2->AppendToTrajectory({q2}, {p2}, {0 * Time::SIUnit()});
    system_.reset(new NBodySystem(
        new std::vector<Body<InertialFrame>*>({body1, body2})));
  }

  SPRKIntegrator integrator_;
  std::unique_ptr<NBodySystem> system_;
};

TEST_F(NBodySystemTest, T) {
  system_->Integrate(integrator_, 10 * Time::SIUnit(), 0.1 * Time::SIUnit(), 1);
}

}  // namespace physics
}  // namespace principia