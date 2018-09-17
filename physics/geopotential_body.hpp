#pragma once

#include "physics/geopotential.hpp"

#include <cmath>

#include "numerics/legendre.hpp"
#include "numerics/polynomial_evaluators.hpp"
#include "quantities/elementary_functions.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace physics {
namespace internal_geopotential {

using numerics::HornerEvaluator;
using numerics::LegendreNormalizationFactor;
using numerics::LegendrePolynomial;
using geometry::InnerProduct;
using quantities::Cos;
using quantities::Inverse;
using quantities::Length;
using quantities::Pow;
using quantities::Sqrt;
using quantities::Square;
using quantities::Sin;
using quantities::SIUnit;

// The notation in this file follows documentation/Geopotential.pdf.
template<typename Frame>
struct Geopotential<Frame>::Precomputations {
  // These quantities are independent from n and m.
  UnitVector x̂;
  UnitVector ŷ;
  UnitVector ẑ;

  Length x;
  Length y;
  Length z;

  Square<Length> r²;
  Length r_norm;

  Angle λ;
  double cos_λ;
  double sin_λ;

  double sin_β;
  double cos_β;

  Vector<Inverse<Length>, Frame> grad_𝔅_vector;
  Vector<Inverse<Length>, Frame> grad_𝔏_vector;

  // These quantities depend on n but are independent from m.
  Inverse<Length> ℜ;
  Vector<Exponentiation<Length, -2>, Frame> grad_ℜ;
};

template<typename Frame>
template<int degree, int order>
struct Geopotential<Frame>::DegreeNOrderM {
  static Vector<Quotient<Acceleration, GravitationalParameter>, Frame>
  Acceleration(OblateBody<Frame> const& body,
               Displacement<Frame> const& r,
               Precomputations const& precomputations);
};

template<typename Frame>
template<int degree, int... orders>
struct Geopotential<Frame>::
DegreeNAllOrders<degree, std::integer_sequence<int, orders...>> {
  static Vector<Quotient<Acceleration, GravitationalParameter>, Frame>
  Acceleration(OblateBody<Frame> const& body,
               Displacement<Frame> const& r,
               Precomputations& precomputations);
};

template<typename Frame>
template<int... degrees>
struct Geopotential<Frame>::AllDegrees<std::integer_sequence<int, degrees...>> {
  static Vector<Quotient<Acceleration, GravitationalParameter>, Frame>
  Acceleration(OblateBody<Frame> const& body,
               Instant const& t,
               Displacement<Frame> const& r,
               Square<Length> const& r²,
               Exponentiation<Length, -3> const& one_over_r³);
};

template<int degree, int order>
double LegendrePolynomialDerivative(double const argument) {
  if constexpr (order > degree) {
    return 0;
  } else {
    static auto const Pn = LegendrePolynomial<degree, HornerEvaluator>();
    static auto const dmpn = Pn.Derivative<order>();
    return dmpn.Evaluate(argument);
  }
}

template<typename Frame>
template<int degree, int order>
Vector<Quotient<Acceleration, GravitationalParameter>, Frame>
Geopotential<Frame>::DegreeNOrderM<degree, order>::Acceleration(
    OblateBody<Frame> const& body,
    Displacement<Frame> const& r,
    Precomputations const& precomputations) {
  if constexpr (degree == 2 && order == 1) {
    return {};
  } else {
    constexpr int n = degree;
    constexpr int m = order;
    static_assert(0 <= m && m <= n);
    static double const normalization_factor =
        LegendreNormalizationFactor(n, m);

    auto const& x̂ = precomputations.x̂;
    auto const& ŷ = precomputations.ŷ;
    auto const& ẑ = precomputations.ẑ;

    auto const& x = precomputations.x;
    auto const& y = precomputations.y;
    auto const& z = precomputations.z;

    auto const& r² = precomputations.r²;
    auto const& r_norm = precomputations.r_norm;

    auto const& λ = precomputations.λ;
    auto const& cos_λ = precomputations.cos_λ;
    auto const& sin_λ = precomputations.sin_λ;

    auto const& cos_β = precomputations.cos_β;
    auto const& sin_β = precomputations.sin_β;

    auto const& grad_𝔅_vector = precomputations.grad_𝔅_vector;
    auto const& grad_𝔏_vector = precomputations.grad_𝔏_vector;

    auto const& ℜ = precomputations.ℜ;
    auto const& grad_ℜ = precomputations.grad_ℜ;

#pragma warning(push)
#pragma warning(disable: 4101)
    double cos_β_to_the_m_minus_1th;  // Not used if m = 0.
#pragma warning(pop)
    double const cos_β_to_the_mth = Pow<m>(cos_β);
    double const Pnm_of_sin_β = LegendrePolynomialDerivative<n, m>(sin_β);
    double const 𝔅 = cos_β_to_the_mth * Pnm_of_sin_β;

    double grad_𝔅_polynomials = cos_β * cos_β_to_the_mth *
                                LegendrePolynomialDerivative<n, m + 1>(sin_β);
    if constexpr (m > 0) {
      cos_β_to_the_m_minus_1th = Pow<m - 1>(cos_β);
      // Remove a singularity when m == 0 and cos_β == 0.
      grad_𝔅_polynomials -= m * sin_β * cos_β_to_the_m_minus_1th * Pnm_of_sin_β;
    }
    Vector<Inverse<Length>, Frame> const grad_𝔅 =
        grad_𝔅_polynomials * grad_𝔅_vector;

    double const Cnm = body.cos()[n][m];
    double const Snm = body.sin()[n][m];

    Angle const mλ = m * λ;
    double const sin_mλ = Sin(mλ);
    double const cos_mλ = Cos(mλ);
    double const 𝔏 = Cnm * cos_mλ + Snm * sin_mλ;

    Vector<Inverse<Length>, Frame> 𝔅_times_grad_𝔏;
    if constexpr (m > 0) {
      // This is not exactly grad_𝔏: we omit the cos_β numerator to remove a
      // singularity.
      Vector<Inverse<Length>, Frame> const grad_𝔏 =
          m * (Snm * cos_mλ - Cnm * sin_mλ) * grad_𝔏_vector;
      // Compensate a cos_β to remove a singularity when cos_β == 0.
      𝔅_times_grad_𝔏 += cos_β_to_the_m_minus_1th * Pnm_of_sin_β * grad_𝔏;
    }

    return normalization_factor *
           (grad_ℜ * 𝔅 * 𝔏 + ℜ * grad_𝔅 * 𝔏 + ℜ * 𝔅_times_grad_𝔏);
  }
}

template<typename Frame>
template<int degree, int... orders>
Vector<Quotient<Acceleration, GravitationalParameter>, Frame> Geopotential<
    Frame>::DegreeNAllOrders<degree, std::integer_sequence<int, orders...>>::
    Acceleration(OblateBody<Frame> const& body,
                 Displacement<Frame> const& r,
                 Precomputations& precomputations) {
  if constexpr (degree < 2) {
    return {};
  } else {
    constexpr int n = degree;

    auto const& r² = precomputations.r²;
    auto const& r_norm = precomputations.r_norm;

    precomputations.ℜ = Pow<n>(body.reference_radius() / r_norm) / r_norm;
    precomputations.grad_ℜ = -(n + 1) * r * precomputations.ℜ / r²;
    return (
        DegreeNOrderM<degree, orders>::Acceleration(body, r, precomputations) +
        ...);
  }
}

template<typename Frame>
template<int... degrees>
Vector<Quotient<Acceleration, GravitationalParameter>, Frame>
Geopotential<Frame>::AllDegrees<std::integer_sequence<int, degrees...>>::
    Acceleration(OblateBody<Frame> const& body,
                 Instant const& t,
                 Displacement<Frame> const& r,
                 Square<Length> const& r²,
                 Exponentiation<Length, -3> const& one_over_r³) {
  auto const from_surface_frame = body.FromSurfaceFrame<SurfaceFrame>(t);
  Precomputations precomputations;

  auto& x̂ = precomputations.x̂;
  auto& ŷ = precomputations.ŷ;
  auto& ẑ = precomputations.ẑ;

  auto& x = precomputations.x;
  auto& y = precomputations.y;
  auto& z = precomputations.z;

  auto& r_norm = precomputations.r_norm;

  auto& λ = precomputations.λ;
  auto& cos_λ = precomputations.cos_λ;
  auto& sin_λ = precomputations.sin_λ;

  auto& cos_β = precomputations.cos_β;
  auto& sin_β = precomputations.sin_β;

  auto& grad_𝔅_vector = precomputations.grad_𝔅_vector;
  auto& grad_𝔏_vector = precomputations.grad_𝔏_vector;

  auto& ℜ = precomputations.ℜ;
  auto& grad_ℜ = precomputations.grad_ℜ;

  x̂ = from_surface_frame(x_);
  ŷ = from_surface_frame(y_);
  ẑ = body.polar_axis();

  x = InnerProduct(r, x̂);
  y = InnerProduct(r, ŷ);
  z = InnerProduct(r, ẑ);

  precomputations.r² = r²;
  r_norm = Sqrt(r²);

  λ = SIUnit<Angle>() * std::atan2(y / SIUnit<Length>(),
                                   x / SIUnit<Length>());
  cos_λ = Cos(λ);
  sin_λ = Sin(λ);

  Square<Length> const x²_plus_y² = x * x + y * y;
  sin_β = z / r_norm;
  cos_β = Sqrt(x²_plus_y²) / r_norm;

  grad_𝔅_vector = (-sin_β * cos_λ * x̂ - sin_β * sin_λ * ŷ + cos_β * ẑ) / r_norm;
  grad_𝔏_vector = (-sin_λ * x̂ + cos_λ * ŷ) / r_norm;

  return (
      DegreeNAllOrders<degrees, std::make_integer_sequence<int, degrees + 1>>::
          Acceleration(body, r, precomputations) +
      ...);
}

template<typename Frame>
Geopotential<Frame>::Geopotential(not_null<OblateBody<Frame> const*> body)
    : body_(body) {}

template<typename Frame>
Vector<Quotient<Acceleration, GravitationalParameter>, Frame>
Geopotential<Frame>::SphericalHarmonicsAcceleration(
    Instant const& t,
    Displacement<Frame> const& r,
    Square<Length> const& r²,
    Exponentiation<Length, -3> const& one_over_r³) const {
  Exponentiation<Length, -2> const one_over_r² = 1 / r²;
  UnitVector const& axis = body_->polar_axis();
  Vector<Quotient<Acceleration, GravitationalParameter>, Frame> acceleration =
      Degree2ZonalAcceleration(axis, r, one_over_r², one_over_r³);
  if (body_->has_c22() || body_->has_s22()) {
    auto const from_surface_frame =
        body_->template FromSurfaceFrame<SurfaceFrame>(t);
    UnitVector const reference = from_surface_frame(x_);
    UnitVector const bireference = from_surface_frame(y_);
    acceleration +=
        Degree2SectoralAcceleration(
            reference, bireference, r, one_over_r², one_over_r³);
  }
  if (body_->has_j3()) {
    acceleration +=
        Degree3ZonalAcceleration(axis, r, r², one_over_r², one_over_r³);
  }
  return acceleration;
}

#define PRINCIPIA_CASE_SPHERICAL_HARMONICS(d)                                  \
  case (d):                                                                    \
    return AllDegrees<std::make_integer_sequence<int, (d + 1)>>::Acceleration( \
        *body_, t, r, r², one_over_r³)

template<typename Frame>
Vector<Quotient<Acceleration, GravitationalParameter>, Frame>
Geopotential<Frame>::GeneralSphericalHarmonicsAcceleration(
    Instant const& t,
    Displacement<Frame> const& r,
    Square<Length> const& r²,
    Exponentiation<Length, -3> const& one_over_r³) const {
  switch (body_->geopotential_degree()) {
    PRINCIPIA_CASE_SPHERICAL_HARMONICS(2);
    PRINCIPIA_CASE_SPHERICAL_HARMONICS(3);
    PRINCIPIA_CASE_SPHERICAL_HARMONICS(4);
    PRINCIPIA_CASE_SPHERICAL_HARMONICS(5);
    case 0:
      return Vector<Quotient<Acceleration, GravitationalParameter>, Frame>{};
    default:
      LOG(FATAL) << "Unexpected degree " << body_->geopotential_degree() << " "
                 << body_->name();
      base::noreturn();
  }
}

#undef PRINCIPIA_CASE_SPHERICAL_HARMONICS

template<typename Frame>
Vector<Quotient<Acceleration, GravitationalParameter>, Frame>
Geopotential<Frame>::Degree2ZonalAcceleration(
    UnitVector const& axis,
    Displacement<Frame> const& r,
    Exponentiation<Length, -2> const& one_over_r²,
    Exponentiation<Length, -3> const& one_over_r³) const {
  Length const r_axis_projection = InnerProduct(axis, r);
  auto const j2_over_r⁵ = body_->j2_over_μ() * one_over_r³ * one_over_r²;
  Vector<Quotient<Acceleration, GravitationalParameter>, Frame> const
      axis_effect = -3 * j2_over_r⁵ * r_axis_projection * axis;
  Vector<Quotient<Acceleration, GravitationalParameter>, Frame> const
      radial_effect =
          j2_over_r⁵ *
          (-1.5 + 7.5 * r_axis_projection * r_axis_projection * one_over_r²) *
          r;
  return axis_effect + radial_effect;
}

template<typename Frame>
Vector<Quotient<Acceleration, GravitationalParameter>, Frame>
Geopotential<Frame>::Degree2SectoralAcceleration(
    UnitVector const& reference,
    UnitVector const& bireference,
    Displacement<Frame> const& r,
    Exponentiation<Length, -2> const& one_over_r²,
    Exponentiation<Length, -3> const& one_over_r³) const {
  Length const r_reference_projection = InnerProduct(reference, r);
  Length const r_bireference_projection = InnerProduct(bireference, r);
  auto const c22_over_r⁵ = body_->c22_over_μ() * one_over_r³ * one_over_r²;
  auto const s22_over_r⁵ = body_->s22_over_μ() * one_over_r³ * one_over_r²;
  Vector<Quotient<Acceleration, GravitationalParameter>, Frame> const
      c22_effect = 6 * c22_over_r⁵ *
                   (-r_bireference_projection * bireference +
                    r_reference_projection * reference +
                    2.5 *
                        (r_bireference_projection * r_bireference_projection -
                         r_reference_projection * r_reference_projection) *
                        one_over_r² * r);
  Vector<Quotient<Acceleration, GravitationalParameter>, Frame> const
      s22_effect = 6 * s22_over_r⁵ *
                   (r_reference_projection * bireference +
                    r_bireference_projection * reference -
                    5 * r_reference_projection * r_bireference_projection *
                        one_over_r² * r);
  return c22_effect + s22_effect;
}

template<typename Frame>
Vector<Quotient<Acceleration, GravitationalParameter>, Frame>
Geopotential<Frame>::Degree3ZonalAcceleration(
    UnitVector const& axis,
    Displacement<Frame> const& r,
    Square<Length> const& r²,
    Exponentiation<Length, -2> const& one_over_r²,
    Exponentiation<Length, -3> const& one_over_r³) const {
  // TODO(phl): Factor the projections across accelerations?
  Length const r_axis_projection = InnerProduct(axis, r);
  Square<Length> const r_axis_projection² =
      r_axis_projection * r_axis_projection;
  auto const j3_over_r⁷ =
      body_->j3_over_μ() * one_over_r³ * one_over_r²* one_over_r²;
  Vector<Quotient<Acceleration, GravitationalParameter>, Frame> const
      axis_effect = 1.5 * j3_over_r⁷ * (r² - 5 * r_axis_projection²) * axis;
  Vector<Quotient<Acceleration, GravitationalParameter>, Frame> const
      radial_effect = j3_over_r⁷ * r_axis_projection *
                      (-7.5 + 17.5 * r_axis_projection² * one_over_r²) * r;
  return axis_effect + radial_effect;
}

template<typename Frame>
const Vector<double, typename Geopotential<Frame>::SurfaceFrame>
    Geopotential<Frame>::x_({1, 0, 0});
template<typename Frame>
const Vector<double, typename Geopotential<Frame>::SurfaceFrame>
    Geopotential<Frame>::y_({0, 1, 0});

}  // namespace internal_geopotential
}  // namespace physics
}  // namespace principia
