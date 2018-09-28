#pragma once

#include "physics/geopotential.hpp"

#include <cmath>

#include "numerics/fixed_arrays.hpp"
#include "numerics/legendre.hpp"
#include "numerics/polynomial_evaluators.hpp"
#include "quantities/elementary_functions.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace physics {
namespace internal_geopotential {

using numerics::FixedVector;
using numerics::HornerEvaluator;
using numerics::LegendreNormalizationFactor;
using numerics::LegendrePolynomial;
using geometry::InnerProduct;
using quantities::ArcTan;
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
template<int size>
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

  double sin_β;
  double cos_β;

  Vector<Inverse<Length>, Frame> grad_𝔅_vector;
  Vector<Inverse<Length>, Frame> grad_𝔏_vector;

  // These quantities depend on n but are independent from m.
  FixedVector<Inverse<Length>, size> ℜ;  // 0 unused.
  Vector<Exponentiation<Length, -2>, Frame> grad_ℜ;

  // These quantities depend on m but are independent from n.
  FixedVector<double, size> cos_mλ;
  FixedVector<double, size> sin_mλ;
  FixedVector<double, size> cos_β_to_the_m;
};

template<typename Frame>
template<int size, int degree, int order>
struct Geopotential<Frame>::DegreeNOrderM {
  static Vector<Quotient<Acceleration, GravitationalParameter>, Frame>
  Acceleration(OblateBody<Frame> const& body,
               Displacement<Frame> const& r,
               Precomputations<size>& precomputations);
};

template<typename Frame>
template<int size, int degree, int... orders>
struct Geopotential<Frame>::
DegreeNAllOrders<size, degree, std::integer_sequence<int, orders...>> {
  static Vector<Quotient<Acceleration, GravitationalParameter>, Frame>
  Acceleration(OblateBody<Frame> const& body,
               Displacement<Frame> const& r,
               Precomputations<size>& precomputations);
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
template<int size, int degree, int order>
Vector<Quotient<Acceleration, GravitationalParameter>, Frame>
Geopotential<Frame>::DegreeNOrderM<size, degree, order>::Acceleration(
    OblateBody<Frame> const& body,
    Displacement<Frame> const& r,
    Precomputations<size>& precomputations) {
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

    auto const& cos_β = precomputations.cos_β;
    auto const& sin_β = precomputations.sin_β;

    auto const& grad_𝔅_vector = precomputations.grad_𝔅_vector;
    auto const& grad_𝔏_vector = precomputations.grad_𝔏_vector;

    auto const& ℜ = precomputations.ℜ[n];
    auto const& grad_ℜ = precomputations.grad_ℜ;

    auto& cos_mλ = precomputations.cos_mλ[m];
    auto& sin_mλ = precomputations.sin_mλ[m];

    auto& cos_β_to_the_m = precomputations.cos_β_to_the_m[m];

    // The fold expressions in the caller ensures that we process n and m by
    // increasing values.  Thus, only the last value of m needs to be
    // initialized for a given value of n.
    if constexpr (m == n) {
      static_assert(m >= 2);

      // Compute the values for m * λ based on the values around m/2 * λ to
      // reduce error accumulation.
      if constexpr (m % 2 == 0) {
        int const h = m / 2;
        double const cos_hλ = precomputations.cos_mλ[h];
        double const sin_hλ = precomputations.sin_mλ[h];
        double const cos_β_to_the_h = precomputations.cos_β_to_the_m[h];
        sin_mλ = 2 * sin_hλ * cos_hλ;
        cos_mλ = cos_hλ * cos_hλ - sin_hλ * sin_hλ;
        cos_β_to_the_m = cos_β_to_the_h * cos_β_to_the_h;
      } else {
        int const h1 = m / 2;
        int const h2 = m - h1;
        double const cos_h1λ = precomputations.cos_mλ[h1];
        double const sin_h1λ = precomputations.sin_mλ[h1];
        double const cos_β_to_the_h1 = precomputations.cos_β_to_the_m[h1];
        double const cos_h2λ = precomputations.cos_mλ[h2];
        double const sin_h2λ = precomputations.sin_mλ[h2];
        double const cos_β_to_the_h2 = precomputations.cos_β_to_the_m[h2];
        sin_mλ = sin_h1λ * cos_h2λ + cos_h1λ * sin_h2λ;
        cos_mλ = cos_h1λ * cos_h2λ - sin_h1λ * sin_h2λ;
        cos_β_to_the_m = cos_β_to_the_h1 * cos_β_to_the_h2;
      }
    }

#pragma warning(push)
#pragma warning(disable: 4101)
    double cos_β_to_the_m_minus_1;  // Not used if m = 0.
#pragma warning(pop)
    double const Pnm_of_sin_β = LegendrePolynomialDerivative<n, m>(sin_β);
    double const 𝔅 = cos_β_to_the_m * Pnm_of_sin_β;

    double grad_𝔅_polynomials = cos_β * cos_β_to_the_m *
                                LegendrePolynomialDerivative<n, m + 1>(sin_β);
    if constexpr (m > 0) {
      cos_β_to_the_m_minus_1 = precomputations.cos_β_to_the_m[m - 1];
      // Remove a singularity when m == 0 and cos_β == 0.
      grad_𝔅_polynomials -= m * sin_β * cos_β_to_the_m_minus_1 * Pnm_of_sin_β;
    }
    Vector<Inverse<Length>, Frame> const grad_𝔅 =
        grad_𝔅_polynomials * grad_𝔅_vector;

    double const Cnm = body.cos()[n][m];
    double const Snm = body.sin()[n][m];
    double const 𝔏 = Cnm * cos_mλ + Snm * sin_mλ;

    Vector<Inverse<Length>, Frame> 𝔅_grad_𝔏;
    if constexpr (m > 0) {
      // This is not exactly grad_𝔏: we omit the cos_β numerator to remove a
      // singularity.
      Vector<Inverse<Length>, Frame> const grad_𝔏 =
          m * (Snm * cos_mλ - Cnm * sin_mλ) * grad_𝔏_vector;
      // Compensate a cos_β to remove a singularity when cos_β == 0.
      𝔅_grad_𝔏 += cos_β_to_the_m_minus_1 * Pnm_of_sin_β * grad_𝔏;
    }

    return normalization_factor *
           (grad_ℜ * 𝔅 * 𝔏 + ℜ * grad_𝔅 * 𝔏 + ℜ * 𝔅_grad_𝔏);
  }
}

template<typename Frame>
template<int size, int degree, int... orders>
Vector<Quotient<Acceleration, GravitationalParameter>, Frame>
Geopotential<Frame>::
DegreeNAllOrders<size, degree, std::integer_sequence<int, orders...>>::
Acceleration(OblateBody<Frame> const& body,
              Displacement<Frame> const& r,
              Precomputations<size>& precomputations) {
  if constexpr (degree < 2) {
    return {};
  } else {
    constexpr int n = degree;

    auto const& r² = precomputations.r²;
    auto const& r_norm = precomputations.r_norm;

    auto& ℜ = precomputations.ℜ[n];
    auto& grad_ℜ = precomputations.grad_ℜ;

    // The fold expressions in the caller ensures that we process n by
    // increasing values.  Thus, we can safely compute ℜ based on values for
    // lower n's.
    if constexpr (n % 2 == 0) {
      int const h = n / 2;
      auto const& ℜh = precomputations.ℜ[h];
      ℜ = ℜh * ℜh * r_norm;
    } else {
      int const h1 = n / 2;
      int const h2 = n - h1;
      auto const& ℜh1 = precomputations.ℜ[h1];
      auto const& ℜh2 = precomputations.ℜ[h2];
      ℜ = ℜh1 * ℜh2 * r_norm;
    }
    grad_ℜ = -(n + 1) * r * ℜ / r²;

    return (... +
            DegreeNOrderM<size, degree, orders>::Acceleration(
                body, r, precomputations));
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
  constexpr int size = sizeof...(degrees);
  auto const from_surface_frame = body.FromSurfaceFrame<SurfaceFrame>(t);
  Precomputations<size> precomputations;

  auto& x̂ = precomputations.x̂;
  auto& ŷ = precomputations.ŷ;
  auto& ẑ = precomputations.ẑ;

  auto& x = precomputations.x;
  auto& y = precomputations.y;
  auto& z = precomputations.z;

  auto& r_norm = precomputations.r_norm;

  auto& cos_β = precomputations.cos_β;
  auto& sin_β = precomputations.sin_β;

  auto& grad_𝔅_vector = precomputations.grad_𝔅_vector;
  auto& grad_𝔏_vector = precomputations.grad_𝔏_vector;

  auto& ℜ1 = precomputations.ℜ[1];

  auto& cos_0λ = precomputations.cos_mλ[0];
  auto& sin_0λ = precomputations.sin_mλ[0];
  auto& cos_1λ = precomputations.cos_mλ[1];
  auto& sin_1λ = precomputations.sin_mλ[1];

  auto& cos_β_to_the_0 = precomputations.cos_β_to_the_m[0];
  auto& cos_β_to_the_1 = precomputations.cos_β_to_the_m[1];

  x̂ = from_surface_frame(x_);
  ŷ = from_surface_frame(y_);
  ẑ = body.polar_axis();

  x = InnerProduct(r, x̂);
  y = InnerProduct(r, ŷ);
  z = InnerProduct(r, ẑ);

  precomputations.r² = r²;
  r_norm = Sqrt(r²);

  Square<Length> const x²_plus_y² = x * x + y * y;
  Length const r_equatorial = Sqrt(x²_plus_y²);

  double cos_λ = 1;
  double sin_λ = 0;
  if (r_equatorial > Length{}) {
    cos_λ = x / r_equatorial;
    sin_λ = y / r_equatorial;
  }

  sin_β = z / r_norm;
  cos_β = r_equatorial / r_norm;

  grad_𝔅_vector = (-sin_β * cos_λ * x̂ - sin_β * sin_λ * ŷ + cos_β * ẑ) / r_norm;
  grad_𝔏_vector = (-sin_λ * x̂ + cos_λ * ŷ) / r_norm;

  ℜ1 = body.reference_radius() / r²;

  cos_0λ = 1;
  sin_0λ = 0;
  cos_1λ = cos_λ;
  sin_1λ = sin_λ;

  cos_β_to_the_0 = 1;
  cos_β_to_the_1 = cos_β;

  // NOTE(phl): The fold expression below should call DegreeNAllOrders with
  // increasing values of the degree.  Unfortunately in VS2017 15.8 it doesn't,
  // in ways that mysteriously depend on the presence of the size parameter in
  // template Precomputations.  We force the right ordering by using an
  // initializer list.
  std::array<Vector<Quotient<quantities::Acceleration,
                             GravitationalParameter>, Frame>,
             size> const accelerations = {
      DegreeNAllOrders<size,
                       degrees,
                       std::make_integer_sequence<int, degrees + 1>>::
          Acceleration(body, r, precomputations)...};

  return (... + accelerations[degrees]);
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
    PRINCIPIA_CASE_SPHERICAL_HARMONICS(6);
    PRINCIPIA_CASE_SPHERICAL_HARMONICS(7);
    PRINCIPIA_CASE_SPHERICAL_HARMONICS(8);
    PRINCIPIA_CASE_SPHERICAL_HARMONICS(9);
    PRINCIPIA_CASE_SPHERICAL_HARMONICS(10);
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
