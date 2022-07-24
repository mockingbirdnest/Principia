
#pragma once

#include "numerics/gradient_descent.hpp"

#include "geometry/grassmann.hpp"
#include "geometry/r3x3_matrix.hpp"
#include "geometry/symmetric_bilinear_form.hpp"
#include "mathematica/mathematica.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace numerics {
namespace internal_gradient_descent {

using geometry::Displacement;
using geometry::InnerProduct;
using geometry::Normalize;
using geometry::R3x3Matrix;
using geometry::SymmetricBilinearForm;
using geometry::SymmetricProduct;
using quantities::Quotient;
using quantities::Square;
using quantities::si::Kilo;
using quantities::si::Metre;
namespace si = quantities::si;

// The type of Dₖ, which approximates the inverse of the Hessian.
template<typename Scalar, typename Frame>
using InverseHessian =
    SymmetricBilinearForm<Quotient<Square<Length>, Scalar>, Frame, Vector>;

constexpr Length initial_step_size = 1 * Metre;//Kilo(Metre);

// Parameters for the Armijo rule.
constexpr double s = 1;
constexpr double β = 0.5;
constexpr double σ = 1e-3;

// We use the BFGS method.
constexpr double ξ = 1;

template<typename Scalar, typename Frame>
double ArmijoRule(Position<Frame> const& xₖ,
                  Displacement<Frame> const& dₖ,
                  Gradient<Scalar, Frame> const& grad_f_xₖ,
                  Field<Scalar, Frame> const& f) {
  Scalar const f_xₖ = f(xₖ);
  Scalar const threshold =  -σ * InnerProduct(grad_f_xₖ, dₖ);
  CHECK_LT(Scalar{}, threshold)
      << "Encountered a Hessian that is not positive definite "
      << "xₖ = " << xₖ << ", dₖ = " << dₖ << " , grad_f_xₖ = " << grad_f_xₖ;
  for (double βᵐs = s;; βᵐs *= β) {
    if (f_xₖ - f(xₖ + βᵐs * dₖ) >= βᵐs * threshold) {
      LOG(ERROR)<<βᵐs;
      return βᵐs;
    }
  }
}

template<typename Scalar, typename Frame>
Position<Frame> GradientDescent(
    Position<Frame> const& start_position,
    Field<Scalar, Frame> const& f,
    Field<Gradient<Scalar, Frame>, Frame> const& grad_f,
    Length const& tolerance) {
  static SymmetricBilinearForm<double, Frame, Vector> const identity =
      SymmetricBilinearForm<double, Frame, Vector>::Identity();

  mathematica::Logger logger(TEMP_DIR / "gradient_descent.wl");

  // The first step uses vanilla steepest descent.
  auto const x₀ = start_position;
  auto const grad_f_x₀ = grad_f(x₀);
  Displacement<Frame> const d₀ = -Normalize(grad_f_x₀) * initial_step_size;
  double const α₀ = ArmijoRule(x₀, d₀, grad_f_x₀, f);
  auto const x₁ = x₀+ α₀ * d₀;

  // Special computation of D₀ using eq. 1.201 and the identity as the starting
  // value of D₀.
  auto const grad_f_x₁ = grad_f(x₁);
  Displacement<Frame> const p₀ = x₁ - x₀;
  auto const q₀ = grad_f_x₁ - grad_f_x₀;
  InverseHessian<Scalar, Frame> const D₀ =
      InnerProduct(p₀, q₀) * identity / q₀.Norm²();

  auto xₖ = x₁;
  auto grad_f_xₖ = grad_f_x₁;
  auto Dₖ = D₀;
  for (;;) {
    LOG(ERROR)<<xₖ;
    LOG(ERROR)<<Dₖ;
    logger.Append("grad",
                  std::tuple{xₖ, grad_f_xₖ},
                  mathematica::ExpressIn(quantities::si::Metre));
    logger.Append("inverseHessian",
                  Dₖ,
                  mathematica::ExpressIn(quantities::si::Metre));
    Displacement<Frame> const dₖ = -Dₖ * grad_f_xₖ;
    LOG(ERROR)<<InnerProduct(grad_f_xₖ, dₖ);
    if (dₖ.Norm() <= tolerance) {
      return xₖ;
    }
    double αₖ = ArmijoRule(xₖ, dₖ, grad_f_xₖ, f);
    Position<Frame> xₖ₊₁;
    Gradient<Scalar, Frame> grad_f_xₖ₊₁;
    do {
      xₖ₊₁ = xₖ + αₖ * dₖ;
      grad_f_xₖ₊₁ = grad_f(xₖ₊₁);
      αₖ *= 1.1;
    LOG(ERROR)<<xₖ₊₁;
      LOG(ERROR)<<-InnerProduct(grad_f_xₖ₊₁, dₖ) <<" "<<
             -0.9 * InnerProduct(grad_f_xₖ, dₖ);
    } while (-InnerProduct(grad_f_xₖ₊₁, dₖ) >
             -0.9 * InnerProduct(grad_f_xₖ, dₖ));
    auto const pₖ = xₖ₊₁ - xₖ;
    auto const qₖ = grad_f_xₖ₊₁ - grad_f_xₖ;
    auto const Dₖqₖ = Dₖ * qₖ;
    auto const τ = InnerProduct(qₖ, Dₖqₖ);
    auto const pₖqₖ = InnerProduct(pₖ, qₖ);
    auto const v = pₖ / pₖqₖ - Dₖqₖ / τ;
    //auto const Dₖ₊₁ = Dₖ + SymmetricProduct(pₖ, pₖ) / pₖqₖ -
    //                  SymmetricProduct(Dₖqₖ, Dₖqₖ) / τ +
    //                  ξ * τ * SymmetricProduct(v, v);
    auto const Dₖ₊₁ =
        Dₖ +
        (pₖqₖ + InnerProduct(qₖ, Dₖqₖ)) * SymmetricProduct(pₖ, pₖ) /
            Pow<2>(pₖqₖ) -
        (SymmetricProduct(Dₖqₖ, pₖ) + SymmetricProduct(pₖ, Dₖqₖ)) / pₖqₖ;

    xₖ = xₖ₊₁;
    grad_f_xₖ = grad_f_xₖ₊₁;
    Dₖ = Dₖ₊₁;
  }
  return xₖ;
}

}  // namespace internal_gradient_descent
}  // namespace numerics
}  // namespace principia
