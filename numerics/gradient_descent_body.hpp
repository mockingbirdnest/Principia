#pragma once

#include "numerics/gradient_descent.hpp"

#include <algorithm>
#include <optional>

#include "absl/status/status.h"
#include "base/stoppable_thread.hpp"  // 🧙 For RETURN_IF_STOPPED.
#include "geometry/symmetric_bilinear_form.hpp"
#include "numerics/elementary_functions.hpp"
#include "numerics/hermite2.hpp"

namespace principia {
namespace numerics {
namespace _gradient_descent {
namespace internal {

using namespace principia::geometry::_symmetric_bilinear_form;
using namespace principia::numerics::_elementary_functions;
using namespace principia::numerics::_hermite2;

// The line search follows [NW06], algorithms 3.5 and 3.6, which guarantee that
// the chosen step obeys the strong Wolfe conditions.

constexpr double c₁ = 1e-4;
constexpr double c₂ = 0.9;
constexpr double α_multiplier = 2;
constexpr double hermite2_tolerance = 0.01;

template<homogeneous_field Scalar, homogeneous_euclidean_space<Scalar> Argument>
  requires(std::totally_ordered<Scalar>)
double Zoom(double α_lo,
            double α_hi,
            Scalar ϕ_α_lo,
            Scalar ϕ_α_hi,
            Scalar ϕʹ_α_lo,
            Scalar const& ϕ_0,
            Scalar const& ϕʹ_0,
            Argument const& x,
            Difference<Argument> const& p,
            VectorField<Scalar, Argument> const& f,
            GateauxDerivative<Scalar, Argument> const& gateaux_derivative_f,
            bool& satisfies_strong_wolfe_condition) {
  std::optional<Scalar> previous_ϕ_αⱼ;
  satisfies_strong_wolfe_condition = true;
  VLOG(1) << "Zoom over: [" << α_lo << " (" << ϕ_α_lo << "), "
                            << α_hi << " (" << ϕ_α_hi << ")]";

#if PRINCIPIA_LOG_ZOOM
  {
    Logger logger(TEMP_DIR / "zoom.wl");
    static constexpr int steps = 100;
    auto const α1 = std::min(α_lo, α_hi);
    auto const α2 = std::max(α_lo, α_hi);
    for (int i = 0; i < steps; ++i) {
      auto const α = α1 + 10 * i * (α2 - α1) / (steps - 1);
      auto const f_α = f(x + α * p);
      logger.Append("phi", std::tuple(α, f_α), ExpressInSIUnits);
    }
    logger.Set("alphaLo", α_lo, ExpressInSIUnits);
    logger.Set("alphaHi", α_hi, ExpressInSIUnits);
    logger.Set("phiAlphaLo", ϕ_α_lo, ExpressInSIUnits);
    logger.Set("phiAlphaHi", ϕ_α_hi, ExpressInSIUnits);
    logger.Set("phiPrimeAlphaLo", ϕʹ_α_lo, ExpressInSIUnits);
    logger.Set("phi0", ϕ_0, ExpressInSIUnits);
    logger.Set("phiPrime0", ϕʹ_0, ExpressInSIUnits);
    logger.Set("x", x, ExpressInSIUnits);
    logger.Set("p", p, ExpressInSIUnits);
  }
#endif

  for (;;) {
    // Note that there is no guarantee here that α_lo < α_hi.
    DCHECK_NE(α_lo, α_hi);
    double αⱼ;
    {
      // Quadratic interpolation.  If the extremum is very close to one of the
      // bounds, zooming would proceed very slowly.  Instead, we bisect, which
      // ensures steady zooming.
      double const α_margin = std::abs(α_hi - α_lo) * hermite2_tolerance;
      Hermite2<Scalar, double> const hermite2(
          {α_lo, α_hi}, {ϕ_α_lo, ϕ_α_hi}, ϕʹ_α_lo);
      auto const α_extremum = hermite2.FindExtremum();
      if (std::min(α_lo, α_hi) + α_margin < α_extremum &&
          α_extremum < std::max(α_lo, α_hi) - α_margin) {
        αⱼ = α_extremum;
      } else {
        // Fall back to bisection.
        αⱼ = (α_lo + α_hi) / 2;
      }
    }

    auto const ϕ_αⱼ = f(x + αⱼ * p);
    VLOG(1) << "  Evaluation at: " << αⱼ << " (" << ϕ_αⱼ << ")";

    // If the function has become (numerically) constant, we might as well
    // return, even though the value of αⱼ may not satisfy the strong Wolfe
    // condition (it probably doesn't, otherwise we would have exited earlier).
    if (previous_ϕ_αⱼ.has_value() && previous_ϕ_αⱼ.value() == ϕ_αⱼ) {
      VLOG(1) << "Numerically constant at: " << αⱼ << " (" << ϕ_αⱼ << ")";
      satisfies_strong_wolfe_condition = false;
      return αⱼ;
    }
    previous_ϕ_αⱼ = ϕ_αⱼ;

    if (ϕ_αⱼ > ϕ_0 + c₁ * αⱼ * ϕʹ_0 || ϕ_αⱼ >= ϕ_α_lo) {
      α_hi = αⱼ;
      ϕ_α_hi = ϕ_αⱼ;
    } else {
      auto const ϕʹ_αⱼ = gateaux_derivative_f(x + αⱼ * p, p);
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

template<homogeneous_field Scalar, homogeneous_euclidean_space<Scalar> Argument>
  requires(std::totally_ordered<Scalar>)
double LineSearch(
    Argument const& x,
    Difference<Argument> const& p,
    Gradient<Scalar, Argument> const& grad_f_x,
    VectorField<Scalar, Argument> const& f,
    GateauxDerivative<Scalar, Argument> const& gateaux_derivative_f,
    bool& satisfies_strong_wolfe_condition) {
  auto const ϕ_0 = f(x);
  auto const ϕʹ_0 = InnerProduct(grad_f_x, p);
  double const αᵢ₋₁ = 0;  // α₀.
  double αᵢ = 1;    // α₁.
  Scalar ϕ_αᵢ₋₁ = ϕ_0;
  Scalar ϕʹ_αᵢ₋₁ = ϕʹ_0;
  satisfies_strong_wolfe_condition = true;
  for (;;) {
    auto const ϕ_αᵢ = f(x + αᵢ * p);
    // For i = 1 the second condition implies the first, so it has no effect on
    // the disjuntion (thus, the test on i in [NW06] is unnecessary).
    if (ϕ_αᵢ > ϕ_0 + c₁ * αᵢ * ϕʹ_0 || ϕ_αᵢ >= ϕ_αᵢ₋₁) {
      return Zoom(αᵢ₋₁, αᵢ,
                  ϕ_αᵢ₋₁, ϕ_αᵢ,
                  ϕʹ_αᵢ₋₁,
                  ϕ_0, ϕʹ_0,
                  x, p, f, gateaux_derivative_f,
                  satisfies_strong_wolfe_condition);
    }
    auto const ϕʹ_αᵢ = gateaux_derivative_f(x + αᵢ * p, p);
    if (Abs(ϕʹ_αᵢ) <= -c₂ * ϕʹ_0) {
      return αᵢ;
    }
    if (ϕʹ_αᵢ >= Scalar{}) {
      return Zoom(αᵢ, αᵢ₋₁,
                  ϕ_αᵢ, ϕ_αᵢ₋₁,
                  ϕʹ_αᵢ,
                  ϕ_0, ϕʹ_0,
                  x, p, f, gateaux_derivative_f,
                  satisfies_strong_wolfe_condition);
    }

    // We don't have a maximum value of α so we blindly increase it until the
    // conditions are met.
    αᵢ *= α_multiplier;
  }
  return αᵢ;
}

template<homogeneous_field Scalar, homogeneous_euclidean_space<Scalar> Argument>
  requires(std::totally_ordered<Scalar>)
absl::StatusOr<Argument> BroydenFletcherGoldfarbShanno(
    Argument const& start_argument,
    VectorField<Scalar, Argument> const& f,
    GradientField<Scalar, Argument> const& grad_f,
    NormType<Difference<Argument>> const& tolerance,
    NormType<Difference<Argument>> const& radius,
    std::optional<NormType<Difference<Argument>>> const&
        first_step) {
  GateauxDerivative<Scalar, Argument> const gateaux_derivative_f =
      [&grad_f](Argument const& argument,
                Difference<Argument> const& direction) {
        return InnerProduct(grad_f(argument), direction);
      };
  return BroydenFletcherGoldfarbShanno(start_argument,
                                       f,
                                       grad_f,
                                       gateaux_derivative_f,
                                       tolerance,
                                       radius,
                                       first_step);
}

// The implementation of BFGS follows [NW06], algorithm 6.18.
template<homogeneous_field Scalar, homogeneous_euclidean_space<Scalar> Argument>
  requires(std::totally_ordered<Scalar>)
absl::StatusOr<Argument> BroydenFletcherGoldfarbShanno(
    Argument const& start_argument,
    VectorField<Scalar, Argument> const& f,
    GradientField<Scalar, Argument> const& grad_f,
    GateauxDerivative<Scalar, Argument> const& gateaux_derivative_f,
    NormType<Difference<Argument>> const& tolerance,
    NormType<Difference<Argument>> const& radius,
    std::optional<NormType<Difference<Argument>>> const&
        first_step) {
  bool satisfies_strong_wolfe_condition;

  // The first step uses vanilla steepest descent.
  auto const x₀ = start_argument;
  auto const grad_f_x₀ = grad_f(x₀);

  VLOG(1) << "Starting from: " << x₀;
  if (grad_f_x₀ == Gradient<Scalar, Argument>{}) {
    VLOG(1) << "Vanishing gradient at: " << x₀;
    return x₀;
  }

  Difference<Argument> const p₀ =
      -Normalize(grad_f_x₀) * first_step.value_or(tolerance);

  double const α₀ = LineSearch(x₀, p₀, grad_f_x₀, f, gateaux_derivative_f,
                               satisfies_strong_wolfe_condition);
  auto const x₁ = x₀+ α₀ * p₀;
  if (!satisfies_strong_wolfe_condition) {
    VLOG(1) << "Doesn't satisfy the strong Wolfe condition at: " << x₁;
    return x₁;
  }

  // Special computation of H₀ using (6.20).
  auto const grad_f_x₁ = grad_f(x₁);
  Difference<Argument> const s₀ = x₁ - x₀;
  auto const y₀ = grad_f_x₁ - grad_f_x₀;
  // TODO(phl): Remove the `internal` use here.
  auto const H₀ = InnerProduct(s₀, y₀) *
                  geometry::_differential_geometry::internal::
                      Generator<Scalar, Argument>::InnerProductForm() /
                  y₀.Norm²();

  auto xₖ = x₁;
  auto grad_f_xₖ = grad_f_x₁;
  auto Hₖ = H₀;
  for (;;) {
    VLOG(1) << "Iterating from: " << xₖ;
    RETURN_IF_STOPPED;
    if ((xₖ - x₀).Norm() > radius) {
      VLOG(1) << "No minimum within search radius at: " << xₖ;
      return absl::Status(termination_condition::NoMinimum, "No minimum found");
    }
    Difference<Argument> const pₖ = -Hₖ * grad_f_xₖ;
    if (pₖ.Norm() <= tolerance) {
      VLOG(1) << "Below tolerance at: " << xₖ
              << ", displacement: " << pₖ.Norm();
      return xₖ;
    }
    double const αₖ = LineSearch(xₖ, pₖ, grad_f_xₖ, f, gateaux_derivative_f,
                                 satisfies_strong_wolfe_condition);
    auto const xₖ₊₁ = xₖ + αₖ * pₖ;
    auto const grad_f_xₖ₊₁ = grad_f(xₖ₊₁);
    auto const sₖ = xₖ₊₁ - xₖ;
    auto const yₖ = grad_f_xₖ₊₁ - grad_f_xₖ;
    auto const sₖyₖ = InnerProduct(sₖ, yₖ);

    // If we can't make progress, e.g., because αₖ is too small, give up.
    if (!satisfies_strong_wolfe_condition) {
      VLOG(1) << "Doesn't satisfy the strong Wolfe condition at: " << xₖ₊₁;
      return xₖ₊₁;
    } else if (sₖyₖ == Scalar{}) {  // NOLINT
      VLOG(1) << "No progress at: " << xₖ₊₁
              << " (s: " << sₖ << ", y: " << yₖ << ")";
      return xₖ₊₁;
    }

    // The formula (6.17) from [NW06] is inconvenient because it uses external
    // products.  Elementary transformations yield the formula below.
    auto const ρ = 1 / sₖyₖ;
    auto const Hₖ₊₁ =
        Hₖ + ρ * ((ρ * Hₖ(yₖ, yₖ) + 1) * SymmetricSquare(sₖ) -
                  2 * SymmetricProduct(Hₖ * yₖ, sₖ));

    xₖ = xₖ₊₁;
    grad_f_xₖ = grad_f_xₖ₊₁;
    Hₖ = Hₖ₊₁;
  }
}

}  // namespace internal
}  // namespace _gradient_descent
}  // namespace numerics
}  // namespace principia
