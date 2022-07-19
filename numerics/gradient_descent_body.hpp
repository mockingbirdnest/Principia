
#pragma once

#include "numerics/gradient_descent.hpp"

#include "geometry/grassmann.hpp"
#include "geometry/symmetric_bilinear_form.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace numerics {
namespace internal_gradient_descent {

using geometry::InnerProduct;
using geometry::SymmetricProduct;
using quantities::si::Kilo;
using quantities::si::Metre;

constexpr Length initial_step_size = 1 * Kilo(Metre);
constexpr double ξ = 1;  // BFGS method.

template<typename Scalar, typename Frame>
Position<Frame> GradientDescent(
    Position<Frame> const& start_position,
    Field<Scalar, Frame> const& f,
    Gradient<Scalar, Frame> const& grad_f,
    TerminationCondition const& termination_condition) {
  Position<Frame> const x₀ = start_position;
  auto const grad_f_x₀ = gradient(x₀);
  Position<Frame> const x₁ = x₀ - initial_step_size * grad_f_x₀;
  auto const grad_f_x₁ = gradient(x₁);
  auto const p₀ = x₁ - x₀;
  auto const q₀ = grad_f_x₁ - grad_f_x₀;
  auto const D₀ = SymmetricProduct(p₀, q₀) / q₀.Norm²();
  while (!termination_condition()) {
    auto const grad_f_xₖ = gradient(xₖ);
    auto const xₖ₊₁ = xₖ - αₖ * Dₖ * grad_f_xₖ;
    auto const pₖ = xₖ₊₁ - xₖ;
    auto const grad_f_xₖ₊₁ = gradient(xₖ₊₁);
    auto const qₖ = grad_f_xₖ₊₁ - grad_f_xₖ;
    auto const Dₖqₖ = Dₖ * qₖ;
    auto const τ = InnerProduct(qₖ, Dₖqₖ);
    auto const pₖqₖ = InnerProduct(pₖ, qₖ);
    auto const v = pₖ / pₖqₖ - Dₖqₖ / τ;
    auto const Dₖ₊₁ = Dₖ + SymmetricProduct(pₖ, pₖ) / pₖqₖ -
                      Dₖ * SymmetricProduct(qₖ, qₖ) * Dₖ / τ +
                      ξ * τ * v.Norm²();
  }
}

}  // namespace internal_gradient_descent
}  // namespace numerics
}  // namespace principia
