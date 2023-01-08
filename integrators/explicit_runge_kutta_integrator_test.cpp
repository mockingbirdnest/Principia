#include "integrators/explicit_runge_kutta_integrator.hpp"

#include <algorithm>
#include <limits>
#include <vector>

#include "base/macros.hpp"
#include "geometry/named_quantities.hpp"
#include "glog/logging.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "quantities/si.hpp"
#include "testing_utilities/almost_equals.hpp"
#include "testing_utilities/approximate_quantity.hpp"
#include "testing_utilities/integration.hpp"
#include "testing_utilities/is_near.hpp"
#include "testing_utilities/matchers.hpp"
#include "testing_utilities/numerics.hpp"

namespace principia {
namespace integrators {
namespace internal_explicit_runge_kutta_integrator {

using geometry::Instant;
using quantities::Abs;
using quantities::Acceleration;
using quantities::AngularFrequency;
using quantities::Cos;
using quantities::Length;
using quantities::Mass;
using quantities::Sin;
using quantities::SpecificImpulse;
using quantities::Speed;
using quantities::Time;
using quantities::si::Centi;
using quantities::si::Kilogram;
using quantities::si::Metre;
using quantities::si::Milli;
using quantities::si::Newton;
using quantities::si::Radian;
using quantities::si::Second;
using testing_utilities::AbsoluteError;
using testing_utilities::AlmostEquals;
using testing_utilities::ComputeHarmonicOscillatorDerivatives1D;
using testing_utilities::EqualsProto;
using testing_utilities::IsNear;
using testing_utilities::StatusIs;
using testing_utilities::operator""_;
using ::std::placeholders::_1;
using ::std::placeholders::_2;
using ::std::placeholders::_3;
using ::testing::ElementsAreArray;
using ::testing::Lt;

using ODE =
    ExplicitFirstOrderOrdinaryDifferentialEquation<Instant, Length, Speed>;

class ExplicitRungeKuttaIntegratorTest
    : public ::testing::Test {};

TEST_F(ExplicitRungeKuttaIntegratorTest,
       Smoke) {
  FixedStepSizeIntegrator<ODE> const& integrator =
      ExplicitRungeKuttaIntegrator<
          methods::RK4, ODE>();
}
}  // namespace internal_ordinary_differential_equations

}  // namespace integrators
}  // namespace principia
