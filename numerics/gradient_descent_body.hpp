#pragma once

#include "numerics/gradient_descent.hpp"

#include <algorithm>
#include <optional>

#include "absl/status/status.h"
#include "base/jthread.hpp"  // 🧙 For RETURN_IF_STOPPED.
#include "geometry/grassmann.hpp"
#include "geometry/point.hpp"
#include "geometry/symmetric_bilinear_form.hpp"
#include "numerics/fixed_arrays.hpp"
#include "numerics/hermite2.hpp"
#include "quantities/elementary_functions.hpp"

namespace principia {
namespace numerics {
namespace _gradient_descent {
namespace internal {

using namespace principia::geometry::_grassmann;
using namespace principia::geometry::_point;
using namespace principia::geometry::_symmetric_bilinear_form;
using namespace principia::numerics::_fixed_arrays;
using namespace principia::numerics::_hermite2;
using namespace principia::quantities::_elementary_functions;

template<typename Scalar, typename S, int s>
struct Generator<Scalar, FixedVector<S, s>> {
  using Gradient = FixedVector<Quotient<Scalar, S>, s>;
  static FixedMatrix<double, s, s> InnerProductForm();
};

template<typename Scalar, typename S, typename F>
struct Generator<Scalar, Vector<S, F>> {
  using Gradient = Vector<Quotient<Scalar, S>, F>;
  static SymmetricBilinearForm<double, F, Vector> InnerProductForm();
};

template<typename Scalar, typename V>
struct Generator<Scalar, Point<V>> {
  using Gradient = typename Generator<Scalar, V>::Gradient;
#if _MSC_FULL_VER == 193'632'532 || \
    _MSC_FULL_VER == 193'632'535 || \
    _MSC_FULL_VER == 193'732'822 || \
    _MSC_FULL_VER == 193'833'135
  using InnerProductFormResult =
      decltype(Generator<Scalar, V>::InnerProductForm());
  static InnerProductFormResult InnerProductForm();
#else
  static decltype(Generator<Scalar, V>::InnerProductForm()) InnerProductForm();
#endif
};

template<typename Scalar, typename S, int s>
FixedMatrix<double, s, s>
Generator<Scalar, FixedVector<S, s>>::InnerProductForm() {
  FixedMatrix<double, s, s> result{};
  for (int i = 0; i < s; ++i) {
    result(i, i) = 1;
  }
  return result;
}

template<typename Scalar, typename S, typename F>
SymmetricBilinearForm<double, F, Vector>
Generator<Scalar, Vector<S, F>>::InnerProductForm() {
  return geometry::_symmetric_bilinear_form::InnerProductForm<F, Vector>();
}

template<typename Scalar, typename V>
#if _MSC_FULL_VER == 193'632'532 || \
    _MSC_FULL_VER == 193'632'535 || \
    _MSC_FULL_VER == 193'732'822 || \
    _MSC_FULL_VER == 193'833'135
typename Generator<Scalar, Point<V>>::InnerProductFormResult
#else
decltype(Generator<Scalar, V>::InnerProductForm())
#endif
Generator<Scalar, Point<V>>::InnerProductForm() {
  return Generator<Scalar, V>::InnerProductForm();
}

// The line search follows [NW06], algorithms 3.5 and 3.6, which guarantee that
// the chosen step obeys the strong Wolfe conditions.

constexpr double c₁ = 1e-4;
constexpr double c₂ = 0.9;
constexpr double α_multiplier = 2;
constexpr double hermite2_tolerance = 0.01;

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
            GateauxDerivative<Scalar, Argument> const& gateaux_derivative_f,
            bool& satisfies_strong_wolfe_condition) {
  std::optional<Scalar> previous_ϕ_αⱼ;
  satisfies_strong_wolfe_condition = true;
  LOG(INFO) << "Zoom over: [" << α_lo << " (" << ϕ_α_lo << "), "
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
    LOG(INFO) << "  Evaluation at: " << αⱼ << " (" << ϕ_αⱼ << ")";

    // If the function has become (numerically) constant, we might as well
    // return, even though the value of αⱼ may not satisfy the strong Wolfe
    // condition (it probably doesn't, otherwise we would have exited earlier).
    if (previous_ϕ_αⱼ.has_value() && previous_ϕ_αⱼ.value() == ϕ_αⱼ) {
      LOG(INFO) << "Numerically constant at: " << αⱼ << " (" << ϕ_αⱼ << ")";
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

template<typename Scalar, typename Argument>
double LineSearch(
    Argument const& x,
    Difference<Argument> const& p,
    Gradient<Scalar, Argument> const& grad_f_x,
    Field<Scalar, Argument> const& f,
    GateauxDerivative<Scalar, Argument> const& gateaux_derivative_f,
    bool& satisfies_strong_wolfe_condition) {
  auto const ϕ_0 = f(x);
  auto const ϕʹ_0 = InnerProduct(grad_f_x, p);
  double αᵢ₋₁ = 0;  // α₀.
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

template<typename Scalar, typename Argument>
absl::StatusOr<Argument> BroydenFletcherGoldfarbShanno(
    Argument const& start_argument,
    Field<Scalar, Argument> const& f,
    Field<Gradient<Scalar, Argument>, Argument> const& grad_f,
    typename Hilbert<Difference<Argument>>::NormType const& tolerance,
    typename Hilbert<Difference<Argument>>::NormType const& radius,
    std::optional<typename Hilbert<Difference<Argument>>::NormType> const&
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
template<typename Scalar, typename Argument>
absl::StatusOr<Argument> BroydenFletcherGoldfarbShanno(
    Argument const& start_argument,
    Field<Scalar, Argument> const& f,
    Field<Gradient<Scalar, Argument>, Argument> const& grad_f,
    GateauxDerivative<Scalar, Argument> const& gateaux_derivative_f,
    typename Hilbert<Difference<Argument>>::NormType const& tolerance,
    typename Hilbert<Difference<Argument>>::NormType const& radius,
    std::optional<typename Hilbert<Difference<Argument>>::NormType> const&
        first_step) {
  bool satisfies_strong_wolfe_condition;

  // The first step uses vanilla steepest descent.
  auto const x₀ = start_argument;
  auto const grad_f_x₀ = grad_f(x₀);

  LOG(INFO) << "Starting from: " << x₀;
  if (grad_f_x₀ == Gradient<Scalar, Argument>{}) {
    LOG(INFO) << "Vanishing gradient at: " << x₀;
    return x₀;
  }

  Difference<Argument> const p₀ =
      -Normalize(grad_f_x₀) * first_step.value_or(tolerance);

  double const α₀ = LineSearch(x₀, p₀, grad_f_x₀, f, gateaux_derivative_f,
                               satisfies_strong_wolfe_condition);
  auto const x₁ = x₀+ α₀ * p₀;
  if (!satisfies_strong_wolfe_condition) {
    LOG(INFO) << "Doesn't satisfy the strong Wolfe condition at: " << x₁;
    return x₁;
  }

  // Special computation of H₀ using (6.20).
  auto const grad_f_x₁ = grad_f(x₁);
  Difference<Argument> const s₀ = x₁ - x₀;
  auto const y₀ = grad_f_x₁ - grad_f_x₀;
  auto const H₀ = InnerProduct(s₀, y₀) *
                  Generator<Scalar, Argument>::InnerProductForm() /
                  y₀.Norm²();

  auto xₖ = x₁;
  auto grad_f_xₖ = grad_f_x₁;
  auto Hₖ = H₀;
  for (;;) {
    LOG(INFO) << "Iterating from: " << xₖ;
    RETURN_IF_STOPPED;
    if ((xₖ - x₀).Norm() > radius) {
      LOG(INFO) << "No minimum within search radius at: " << xₖ;
      return absl::Status(termination_condition::NoMinimum, "No minimum found");
    }
    Difference<Argument> const pₖ = -Hₖ * grad_f_xₖ;
    if (pₖ.Norm() <= tolerance) {
      LOG(INFO) << "Below tolerance at: " << xₖ
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
      LOG(INFO) << "Doesn't satisfy the strong Wolfe condition at: " << xₖ₊₁;
      return xₖ₊₁;
    } else if (sₖyₖ == Scalar{}) {  // NOLINT
      LOG(INFO) << "No progress at: " << xₖ₊₁
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
