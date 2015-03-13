
#include "integrators/symplectic_partitioned_runge_kutta_integrator.hpp"

#include <algorithm>
#include <string>
#include <vector>

#include "geometry/sign.hpp"
#include "glog/logging.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "quantities/quantities.hpp"
#include "quantities/named_quantities.hpp"
#include "testing_utilities/almost_equals.hpp"
#include "testing_utilities/numerical_analysis.hpp"
#include "testing_utilities/numerics.hpp"
#include "testing_utilities/statistics.hpp"

namespace principia {

using quantities::Abs;
using quantities::Acceleration;
using quantities::AngularFrequency;
using quantities::Energy;
using quantities::Force;
using quantities::Length;
using quantities::Mass;
using quantities::Speed;
using quantities::Pow;
using quantities::Power;
using quantities::SIUnit;
using quantities::Speed;
using quantities::Stiffness;
using quantities::Time;
using si::Joule;
using si::Kilogram;
using si::Metre;
using si::Newton;
using si::Second;
using testing_utilities::AlmostEquals;
using testing_utilities::BidimensionalDatasetMathematicaInput;
using testing_utilities::ComputeHarmonicOscillatorAcceleration;
using testing_utilities::PearsonProductMomentCorrelationCoefficient;
using testing_utilities::RelativeError;
using testing_utilities::Slope;
using ::testing::AllOf;
using ::testing::Eq;
using ::testing::Gt;
using ::testing::Lt;
using ::testing::Ne;
using ::testing::ValuesIn;

#define INTEGRATOR(name) &name(), #name

namespace integrators {

namespace {

struct SRKN {
  not_null<SRKNIntegrator const*> integrator;
  std::string name;
};

// This allows the test output to be legible, i.e.,
// "where GetParam() = Leapfrog" rather than
// "where GetParam() = n-byte object <hex>"
std::ostream& operator<<(std::ostream& stream, SRKN param) {
  return stream << param.name;
}

std::vector<SRKN> Instances() {
  return {
      {INTEGRATOR(Leapfrog)},
      {INTEGRATOR(PseudoLeapfrog)},
      {INTEGRATOR(McLachlanAtela1992Order2Optimal)},
      {INTEGRATOR(Ruth1983)},
      {INTEGRATOR(McLachlanAtela1992Order3Optimal)},
      {INTEGRATOR(CandyRozmus1991ForestRuth1990SynchronousMomenta)},
      {INTEGRATOR(CandyRozmus1991ForestRuth1990SynchronousPositions)},
      {INTEGRATOR(McLachlanAtela1992Order4Optimal)},
      {INTEGRATOR(McLachlanAtela1992Order5Optimal)},
      {INTEGRATOR(Yoshida1990Order6A)},
      {INTEGRATOR(Yoshida1990Order6B)},
      {INTEGRATOR(Yoshida1990Order6C)},
      {INTEGRATOR(Yoshida1990Order8A)},
      {INTEGRATOR(Yoshida1990Order8B)},
      {INTEGRATOR(Yoshida1990Order8C)},
      {INTEGRATOR(Yoshida1990Order8D)},
      {INTEGRATOR(Yoshida1990Order8E)}};
}

}  // namespace

class SRKNTest : public testing::TestWithParam<SRKN> {
 public:
  static void SetUpTestCase() {
    google::LogToStderr();
  }

 protected:
  SRKNTest() : integrator_(GetParam().integrator) {}

  not_null<SRKNIntegrator const*> const     integrator_;
  SRKNIntegrator::Parameters<Length, Speed> parameters_;
  SRKNIntegrator::Solution<Length, Speed>   solution_;
};

INSTANTIATE_TEST_CASE_P(SRKNTests, SRKNTest, ValuesIn(Instances()));

TEST_P(SRKNTest, ConsistentWeights) {
  // Check that the time argument of the force computation is correct by
  // integrating uniform linear motion.
  // We check this for all schemes.
  Speed const v = 1 * Metre / Second;
  Mass const m = 1 * Kilogram;
  auto compute_acceleration = [v](
      Time const& t,
      std::vector<Length> const& q,
      not_null<std::vector<Acceleration>*> const result) {
    EXPECT_THAT(q[0], AlmostEquals(v * t, 0, 4096));
    (*result)[0] = 0 * Metre / Pow<2>(Second);
  };
  parameters_.initial.positions.emplace_back(0 * Metre);
  parameters_.initial.momenta.emplace_back(v);
  parameters_.initial.time = 0 * Second;
  parameters_.tmax = 16 * Second;
  parameters_.Δt = 1 * Second;
  parameters_.sampling_period = 5;
  integrator_->SolveTrivialKineticEnergyIncrement<Length, Speed>(
      compute_acceleration, parameters_, &solution_);
  EXPECT_THAT(solution_.back().positions.back().value,
              AlmostEquals(v * parameters_.tmax, 0, 4));
}

TEST_P(SRKNTest, ExactInexactTMax) {
  parameters_.initial.positions.emplace_back(SIUnit<Length>());
  parameters_.initial.momenta.emplace_back(Speed());
  parameters_.initial.time = Time();
  parameters_.tmax = 10.0 * SIUnit<Time>();
  parameters_.sampling_period = 1;
  parameters_.Δt = (1.0 / 3.000001) * SIUnit<Time>();
  parameters_.tmax_is_exact = false;
  integrator_->SolveTrivialKineticEnergyIncrement<Length, Speed>(
      &ComputeHarmonicOscillatorAcceleration,
      parameters_,
      &solution_);
  EXPECT_EQ(30, solution_.size());
  EXPECT_THAT(solution_.back().time.value, Lt(parameters_.tmax));
  EXPECT_THAT(solution_.back().time.error, Ne(0.0 * SIUnit<Time>()));

  parameters_.tmax_is_exact = true;
  integrator_->SolveTrivialKineticEnergyIncrement<Length, Speed>(
      &ComputeHarmonicOscillatorAcceleration,
      parameters_,
      &solution_);
  EXPECT_EQ(30, solution_.size());
  EXPECT_THAT(solution_.back().time.value, Eq(parameters_.tmax));
  EXPECT_THAT(solution_.back().time.error, Eq(0.0 * SIUnit<Time>()));

  parameters_.Δt = (1.0 / 2.999999) * SIUnit<Time>();
  parameters_.tmax_is_exact = false;
  integrator_->SolveTrivialKineticEnergyIncrement<Length, Speed>(
      &ComputeHarmonicOscillatorAcceleration,
      parameters_,
      &solution_);
  EXPECT_EQ(29, solution_.size());
  EXPECT_THAT(solution_.back().time.value, Lt(parameters_.tmax));
  EXPECT_THAT(solution_.back().time.error, Ne(0.0 * SIUnit<Time>()));

  parameters_.tmax_is_exact = true;
  integrator_->SolveTrivialKineticEnergyIncrement<Length, Speed>(
      &ComputeHarmonicOscillatorAcceleration,
      parameters_,
      &solution_);
  EXPECT_EQ(30, solution_.size());
  EXPECT_THAT(solution_.back().time.value, Eq(parameters_.tmax));
  EXPECT_THAT(solution_.back().time.error, Eq(0.0 * SIUnit<Time>()));

  parameters_.Δt = 11.0 * SIUnit<Time>();
  parameters_.tmax_is_exact = false;
  integrator_->SolveTrivialKineticEnergyIncrement<Length, Speed>(
      &ComputeHarmonicOscillatorAcceleration,
      parameters_,
      &solution_);
  EXPECT_EQ(0, solution_.size());

  parameters_.tmax_is_exact = true;
  integrator_->SolveTrivialKineticEnergyIncrement<Length, Speed>(
      &ComputeHarmonicOscillatorAcceleration,
      parameters_,
      &solution_);
  EXPECT_EQ(1, solution_.size());
  EXPECT_THAT(solution_.back().time.value, Eq(parameters_.tmax));
  EXPECT_THAT(solution_.back().time.error, Eq(0.0 * SIUnit<Time>()));

  parameters_.Δt = 100.0 * SIUnit<Time>();
  parameters_.tmax_is_exact = false;
  integrator_->SolveTrivialKineticEnergyIncrement<Length, Speed>(
      &ComputeHarmonicOscillatorAcceleration,
      parameters_,
      &solution_);
  EXPECT_EQ(0, solution_.size());

  parameters_.tmax_is_exact = true;
  integrator_->SolveTrivialKineticEnergyIncrement<Length, Speed>(
      &ComputeHarmonicOscillatorAcceleration,
      parameters_,
      &solution_);
  EXPECT_EQ(1, solution_.size());
  EXPECT_THAT(solution_.back().time.value, Eq(parameters_.tmax));
  EXPECT_THAT(solution_.back().time.error, Eq(0.0 * SIUnit<Time>()));
}

}  // namespace integrators
}  // namespace principia
