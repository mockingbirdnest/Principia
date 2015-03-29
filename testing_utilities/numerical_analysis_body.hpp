#pragma once

#include <vector>

#include "quantities/elementary_functions.hpp"
#include "quantities/quantities.hpp"
#include "quantities/named_quantities.hpp"

namespace principia {

using geometry::Displacement;
using geometry::InnerProduct;
using quantities::Force;
using quantities::Exponentiation;
using quantities::GravitationalParameter;
using quantities::Length;
using quantities::Mass;
using quantities::Momentum;
using quantities::Pow;
using quantities::SIUnit;
using quantities::Speed;
using quantities::Sqrt;
using quantities::Stiffness;
using quantities::Time;

namespace testing_utilities {

inline void ComputeHarmonicOscillatorForce(
    Time const& t,
    std::vector<Length> const& q,
    std::vector<Force>* const result) {
  (*result)[0] = -q[0] * SIUnit<Stiffness>();
}

inline void ComputeHarmonicOscillatorVelocity(
    std::vector<Momentum> const& p,
    std::vector<Speed>* const result) {
  (*result)[0] = p[0] / SIUnit<Mass>();
}

inline void ComputeHarmonicOscillatorAcceleration(
    Time const& t,
    std::vector<Length> const& q,
    std::vector<Acceleration>* const result) {
  (*result)[0] = -q[0] * (SIUnit<Stiffness>() / SIUnit<Mass>());
}

inline void ComputeKeplerAcceleration(
    Time const& t,
    std::vector<Length> const& q,
    std::vector<Acceleration>* const result) {
  auto const r_squared = q[0] * q[0] + q[1] * q[1];
  auto const minus_μ_over_r_cubed =
      -SIUnit<GravitationalParameter>() * Sqrt(r_squared) /
          (r_squared * r_squared);
  (*result)[0] = q[0] * minus_μ_over_r_cubed;
  (*result)[1] = q[1] * minus_μ_over_r_cubed;
}

template<typename Frame>
void ComputeGravitationalAcceleration(
    Time const& t,
    std::vector<Position<Frame>> const& q,
    std::vector<Vector<Acceleration, Frame>>* const result,
    std::vector<MassiveBody> const& bodies) {
  result->assign(result->size(), Vector<Acceleration, Frame>());
  for (int b1 = 1; b1 < q.size(); ++b1) {
    GravitationalParameter const& μ1 = bodies[b1].gravitational_parameter();
    for (int b2 = 0; b2 < b1; ++b2) {
      Displacement<Frame> const Δq = q[b1] - q[b2];
      Exponentiation<Length, 2> const r_squared = InnerProduct(Δq, Δq);
      Exponentiation<Length, -3> const one_over_r_cubed =
          Sqrt(r_squared) / (r_squared * r_squared);
      {
        auto const μ1_over_r_cubed = μ1 * one_over_r_cubed;
        (*result)[b2] += Δq * μ1_over_r_cubed;
      }
      // Lex. III. Actioni contrariam semper & æqualem esse reactionem:
      // sive corporum duorum actiones in se mutuo semper esse æquales &
      // in partes contrarias dirigi.
      {
        GravitationalParameter const& μ2 = bodies[b2].gravitational_parameter();
        auto const μ2_over_r_cubed = μ2 * one_over_r_cubed;
        (*result)[b1] -= Δq * μ2_over_r_cubed;
      }
    }
  }
}

}  // namespace testing_utilities
}  // namespace principia
