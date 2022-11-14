
#pragma once

#include "numerics/gradient_descent.hpp"

#include <optional>

#include "geometry/grassmann.hpp"
#include "geometry/r3x3_matrix.hpp"
#include "geometry/symmetric_bilinear_form.hpp"
#include "numerics/hermite2.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace numerics {
namespace internal_gradient_descent {

using geometry::InnerProduct;
using geometry::InnerProductForm;
using geometry::Normalize;
using geometry::SymmetricBilinearForm;
using geometry::SymmetricProduct;
using geometry::Vector;
using quantities::Abs;
using quantities::Quotient;
using quantities::Square;
using quantities::si::Metre;
namespace si = quantities::si;

// A helper to use |Argument| with |SymmetricBilinearForm|.
template<typename A>
struct ArgumentHelper;

template<typename Scalar, typename Frame>
struct ArgumentHelper<Vector<Scalar, Frame>> {
  static SymmetricBilinearForm<double, Frame, Vector> InnerProductForm() {
    return geometry::InnerProductForm<Frame, Vector>();
  }
};

// The line search follows [NW06], algorithms 3.5 and 3.6, which guarantee that
// the chosen step obeys the strong Wolfe conditions.

constexpr double c₁ = 1e-4;
constexpr double c₂ = 0.9;
constexpr double α_multiplier = 2;

template<typename Scalar, typename Argument>
double Zoom(double α_lo,
            double α_hi,
            Scalar ϕ_α_lo,
            Scalar ϕ_α_hi,
            Scalar ϕʹ_α_lo,
            Scalar const& ϕ_0,
            Scalar const& ϕʹ_0,
            Argument const& x,
            Difference<Argument> const& p,
            Field<Scalar, Argument> const& f,
            Field<Gradient<Scalar, Argument>, Argument> const& grad_f) {
  std::optional<Scalar> previous_ϕ_αⱼ;
  for (;;) {
    DCHECK_NE(α_lo, α_hi);
    double αⱼ;
    {
      // Quadratic interpolation.
      Hermite2<Scalar, double> const hermite2(
          {α_lo, α_hi}, {ϕ_α_lo, ϕ_α_hi}, ϕʹ_α_lo);
      auto const α_extremum = hermite2.FindExtremum();
      if (α_lo < α_extremum && α_extremum < α_hi) {
        αⱼ = α_extremum;
      } else {
        // Fall back to bisection.
        αⱼ = (α_lo + α_hi) / 2;
      }
    }

    auto const ϕ_αⱼ = f(x + αⱼ * p);

    // If the function has become (numerically) constant, we might as well
    // return, even though the value of αⱼ may not meet the strong Wolfe
    // condition.
    if (previous_ϕ_αⱼ.has_value() && previous_ϕ_αⱼ.value() == ϕ_αⱼ) {
      return αⱼ;
    }
    previous_ϕ_αⱼ = ϕ_αⱼ;

    if (ϕ_αⱼ > ϕ_0 + c₁ * αⱼ * ϕʹ_0 || ϕ_αⱼ >= ϕ_α_lo) {
      α_hi = αⱼ;
      ϕ_α_hi = ϕ_αⱼ;
    } else {
      auto const ϕʹ_αⱼ = InnerProduct(grad_f(x + αⱼ * p), p);
      if (Abs(ϕʹ_αⱼ) <= -c₂ * ϕʹ_0) {
        return αⱼ;
      }
      if (ϕʹ_αⱼ * (α_hi - α_lo) >= Scalar{}) {
        α_hi = α_lo;
        ϕ_α_hi = ϕ_α_lo;
      }
      α_lo = αⱼ;
      ϕ_α_lo = ϕ_αⱼ;
      ϕʹ_α_lo = ϕʹ_αⱼ;
    }
  }
}

template<typename Scalar, typename Argument>
double LineSearch(Argument const& x,
                  Difference<Argument> const& p,
                  Field<Scalar, Argument> const& f,
                  Field<Gradient<Scalar, Argument>, Argument> const& grad_f) {
  auto const ϕ_0 = f(x);
  auto const ϕʹ_0 = InnerProduct(grad_f(x), p);
  double αᵢ₋₁ = 0;  // α₀.
  double αᵢ = 1;    // α₁.
  Scalar ϕ_αᵢ₋₁ = ϕ_0;
  Scalar ϕʹ_αᵢ₋₁ = ϕʹ_0;
  for (;;) {
    auto const ϕ_αᵢ = f(x + αᵢ * p);
    // For i = 1 the second condition implies the first, so it has no effect on
    // the disjuntion (thus, the test on i in [NW06] is unnecessary).
    if (ϕ_αᵢ > ϕ_0 + c₁ * αᵢ * ϕʹ_0 || ϕ_αᵢ >= ϕ_αᵢ₋₁) {
      return Zoom(αᵢ₋₁, αᵢ,
                  ϕ_αᵢ₋₁, ϕ_αᵢ,
                  ϕʹ_αᵢ₋₁,
                  ϕ_0, ϕʹ_0,
                  x, p, f, grad_f);
    }
    auto const ϕʹ_αᵢ = InnerProduct(grad_f(x + αᵢ * p), p);
    if (Abs(ϕʹ_αᵢ) <= -c₂ * ϕʹ_0) {
      return αᵢ;
    }
    if (ϕʹ_αᵢ >= Scalar{}) {
      return Zoom(αᵢ, αᵢ₋₁,
                  ϕ_αᵢ, ϕ_αᵢ₋₁,
                  ϕʹ_αᵢ,
                  ϕ_0, ϕʹ_0,
                  x, p, f, grad_f);
    }

    // We don't have a maximum value of α so we blindly increase it until the
    // conditions are met.
    αᵢ *= α_multiplier;
  }
  return αᵢ;
}

// The implementation of BFGS follows [NW06], algorithm 6.18.
template<typename Scalar, typename Argument>
Argument BroydenFletcherGoldfarbShanno(
    Argument const& start_argument,
    Field<Scalar, Argument> const& f,
    Field<Gradient<Scalar, Argument>, Argument> const& grad_f,
    Length const& tolerance) {
  // The first step uses vanilla steepest descent.
  auto const x₀ = start_argument;
  auto const grad_f_x₀ = grad_f(x₀);

  if (grad_f_x₀ == Gradient<Scalar, Argument>{}) {
    return x₀;
  }

  // We (ab)use the tolerance to determine the first step size.  The assumption
  // is that, if the caller provides a reasonable value then (1) we won't miss
  // "interesting features" of f; (2) the finite differences won't underflow or
  // have other unpleasant properties.
  Difference<Argument> const p₀ = -Normalize(grad_f_x₀) * tolerance;

  double const α₀ = LineSearch(x₀, p₀, f, grad_f);
  auto const x₁ = x₀+ α₀ * p₀;

  // Special computation of H₀ using (6.20).
  auto const grad_f_x₁ = grad_f(x₁);
  Difference<Argument> const s₀ = x₁ - x₀;
  auto const y₀ = grad_f_x₁ - grad_f_x₀;
  auto const H₀ = InnerProduct(s₀, y₀) *
                  ArgumentHelper<Difference<Argument>>::InnerProductForm() /
                  y₀.Norm²();

  auto xₖ = x₁;
  auto grad_f_xₖ = grad_f_x₁;
  auto Hₖ = H₀;
  for (;;) {
    Difference<Argument> const pₖ = -Hₖ * grad_f_xₖ;
    if (pₖ.Norm() <= tolerance) {
      return xₖ;
    }
    double const αₖ = LineSearch(xₖ, pₖ, f, grad_f);
    auto const xₖ₊₁ = xₖ + αₖ * pₖ;
    auto const grad_f_xₖ₊₁ = grad_f(xₖ₊₁);

    auto const sₖ = xₖ₊₁ - xₖ;
    auto const yₖ = grad_f_xₖ₊₁ - grad_f_xₖ;
    auto const sₖyₖ = InnerProduct(sₖ, yₖ);

    // If we can't make progress, e.g., because αₖ is too small, give up.
    if (sₖyₖ == Scalar{}) {
      return xₖ;
    }

    // The formula (6.17) from [NW06] is inconvenient because it uses external
    // products.  Elementary transformations yield the formula below.
    auto const ρ = 1 / sₖyₖ;
    auto const Hₖ₊₁ =
        Hₖ + ρ * ((ρ * Hₖ(yₖ, yₖ) + 1) * SymmetricProduct(sₖ, sₖ) -
                  2 * SymmetricProduct(Hₖ * yₖ, sₖ));

    xₖ = xₖ₊₁;
    grad_f_xₖ = grad_f_xₖ₊₁;
    Hₖ = Hₖ₊₁;
  }
  return xₖ;
}

}  // namespace internal_gradient_descent
}  // namespace numerics
}  // namespace principia
