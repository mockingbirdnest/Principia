
#pragma once

#include <vector>

#include "testing_utilities/integration.hpp"

#include "quantities/elementary_functions.hpp"
#include "quantities/quantities.hpp"
#include "quantities/named_quantities.hpp"

namespace principia {
namespace testing_utilities {
namespace internal_integration {

using geometry::Displacement;
using geometry::InnerProduct;
using quantities::Exponentiation;
using quantities::Force;
using quantities::GravitationalParameter;
using quantities::Length;
using quantities::Mass;
using quantities::Momentum;
using quantities::Pow;
using quantities::SIUnit;
using quantities::Speed;
using quantities::Sqrt;
using quantities::Square;
using quantities::Stiffness;
using quantities::Time;

inline void ComputeHarmonicOscillatorForce(
    Time const& t,
    std::vector<Length> const& q,
    std::vector<Force>& result) {
  result[0] = -q[0] * SIUnit<Stiffness>();
}

inline void ComputeHarmonicOscillatorVelocity(
    std::vector<Momentum> const& p,
    std::vector<Speed>& result) {
  result[0] = p[0] / SIUnit<Mass>();
}

inline Status ComputeHarmonicOscillatorAcceleration(
    Instant const& t,
    std::vector<Length> const& q,
    std::vector<Acceleration>& result,
    int* evaluations) {
  result[0] = -q[0] * (SIUnit<Stiffness>() / SIUnit<Mass>());
  if (evaluations != nullptr) {
    ++*evaluations;
  }
  return Status::OK;
}

inline Status ComputeKeplerAcceleration(
    Instant const& t,
    std::vector<Length> const& q,
    std::vector<Acceleration>& result,
    int* evaluations) {
  auto const r² = q[0] * q[0] + q[1] * q[1];
  auto const minus_μ_over_r³ =
      -SIUnit<GravitationalParameter>() * Sqrt(r²) / (r² * r²);
  result[0] = q[0] * minus_μ_over_r³;
  result[1] = q[1] * minus_μ_over_r³;
  if (evaluations != nullptr) {
    ++*evaluations;
  }
  return Status::OK;
}

template<typename Frame>
void ComputeGravitationalAcceleration(
    Time const& t,
    std::vector<Position<Frame>> const& q,
    std::vector<Vector<Acceleration, Frame>>& result,
    std::vector<MassiveBody> const& bodies) {
  result.assign(result.size(), Vector<Acceleration, Frame>());
  for (int b1 = 1; b1 < q.size(); ++b1) {
    GravitationalParameter const& μ1 = bodies[b1].gravitational_parameter();
    for (int b2 = 0; b2 < b1; ++b2) {
      Displacement<Frame> const Δq = q[b1] - q[b2];
      Square<Length> const r² = Δq.Norm²();
      Exponentiation<Length, -3> const one_over_r³ = Sqrt(r²) / (r² * r²);
      {
        auto const μ1_over_r³ = μ1 * one_over_r³;
        result[b2] += Δq * μ1_over_r³;
      }
      // Lex. III. Actioni contrariam semper & æqualem esse reactionem:
      // sive corporum duorum actiones in se mutuo semper esse æquales &
      // in partes contrarias dirigi.
      {
        GravitationalParameter const& μ2 = bodies[b2].gravitational_parameter();
        auto const μ2_over_r³ = μ2 * one_over_r³;
        result[b1] -= Δq * μ2_over_r³;
      }
    }
  }
}

}  // namespace internal_integration
}  // namespace testing_utilities
}  // namespace principia
