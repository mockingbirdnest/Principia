#pragma once

#include "physics/geopotential.hpp"

#include <cmath>

#include "geometry/grassmann.hpp"
#include "geometry/r3_element.hpp"
#include "numerics/fixed_arrays.hpp"
#include "numerics/legendre.hpp"
#include "numerics/polynomial_evaluators.hpp"
#include "quantities/elementary_functions.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace physics {
namespace internal_geopotential {

using numerics::FixedLowerTriangularMatrix;
using numerics::FixedVector;
using numerics::HornerEvaluator;
using numerics::LegendreNormalizationFactor;
using numerics::LegendrePolynomial;
using geometry::InnerProduct;
using geometry::R3Element;
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
  typename OblateBody<Frame>::GeopotentialCoefficients const* cos;
  typename OblateBody<Frame>::GeopotentialCoefficients const* sin;

  Square<Length> r²;
  Length r_norm;

  double sin_β;
  double cos_β;

  Vector<Inverse<Length>, SurfaceFrame> grad_𝔅_vector;
  Vector<Inverse<Length>, SurfaceFrame> grad_𝔏_vector;

  // These quantities depend on n but are independent from m.
  FixedVector<Inverse<Length>, size> ℜ;  // 0 unused.
  Vector<Exponentiation<Length, -2>, SurfaceFrame> grad_ℜ;

  // These quantities depend on m but are independent from n.
  FixedVector<double, size> cos_mλ;
  FixedVector<double, size> sin_mλ;
  FixedVector<double, size> cos_β_to_the_m;

  // These quantities depend on both n and m.  Note that the zeros for m > n are
  // not stored.
  FixedLowerTriangularMatrix<double, size> DmPn_of_sin_β;
};

template<typename Frame>
template<int size, int degree, int order>
struct Geopotential<Frame>::DegreeNOrderM {
  FORCE_INLINE(static)
  auto Acceleration(Precomputations<size>& precomputations)
      -> Vector<ReducedAcceleration, SurfaceFrame>;
};

template<typename Frame>
template<int size, int degree, int... orders>
struct Geopotential<Frame>::
DegreeNAllOrders<size, degree, std::integer_sequence<int, orders...>> {
  static auto Acceleration(Displacement<SurfaceFrame> const& r_surface,
                           Precomputations<size>& precomputations)
      -> Vector<ReducedAcceleration, SurfaceFrame>;
};

template<typename Frame>
template<int... degrees>
struct Geopotential<Frame>::AllDegrees<std::integer_sequence<int, degrees...>> {
  static Vector<ReducedAcceleration, Frame>
  Acceleration(OblateBody<Frame> const& body,
               Instant const& t,
               Displacement<Frame> const& r,
               Square<Length> const& r²,
               Exponentiation<Length, -3> const& one_over_r³);
};

template<typename Frame>
template<int size, int degree, int order>
auto Geopotential<Frame>::DegreeNOrderM<size, degree, order>::Acceleration(
    Precomputations<size>& precomputations)
    -> Vector<ReducedAcceleration, SurfaceFrame> {
  if constexpr (degree == 2 && order == 1) {
    // Let's not forget the Legendre derivative that we would compute if we did
    // not short-circuit.
    precomputations.DmPn_of_sin_β[2][2] = 3;
    return {};
  } else {
    constexpr int n = degree;
    constexpr int m = order;
    static_assert(0 <= m && m <= n);
    static double const normalization_factor =
        LegendreNormalizationFactor(n, m);

    auto const& cos_β = precomputations.cos_β;
    auto const& sin_β = precomputations.sin_β;

    auto const& grad_𝔅_vector = precomputations.grad_𝔅_vector;
    auto const& grad_𝔏_vector = precomputations.grad_𝔏_vector;

    auto const& ℜ = precomputations.ℜ[n];
    auto const& grad_ℜ = precomputations.grad_ℜ;

    auto& cos_mλ = precomputations.cos_mλ[m];
    auto& sin_mλ = precomputations.sin_mλ[m];

    auto& cos_β_to_the_m = precomputations.cos_β_to_the_m[m];

    auto& DmPn_of_sin_β = precomputations.DmPn_of_sin_β;
    auto const& cos = precomputations.cos;
    auto const& sin = precomputations.sin;

    // The caller ensures that we process n and m by increasing values.  Thus,
    // only the last value of m needs to be initialized for a given value of n.
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

    // Recurrence relationship between the Legendre polynomials.
    if constexpr (m == 0) {
      static_assert(n >= 2);
      DmPn_of_sin_β[n][0] = ((2 * n - 1) * sin_β * DmPn_of_sin_β[n - 1][0] -
                             (n - 1) * DmPn_of_sin_β[n - 2][0]) /
                            n;
    }

    // Recurrence relationship between the associated Legendre polynomials.
    // Account for the fact that DmPn_of_sin_β is identically zero if m > n.
    if constexpr (m == n) {
      // Do not store the zero.
    } else if constexpr (m == n - 1) {  // NOLINT(readability/braces)
      static_assert(n >= 1);
      DmPn_of_sin_β[n][m + 1] =
          ((2 * n - 1) * (m + 1) * DmPn_of_sin_β[n - 1][m]) / n;
    } else if constexpr (m == n - 2) {  // NOLINT(readability/braces)
      static_assert(n >= 1);
      DmPn_of_sin_β[n][m + 1] =
          ((2 * n - 1) * (sin_β * DmPn_of_sin_β[n - 1][m + 1] +
                          (m + 1) * DmPn_of_sin_β[n - 1][m])) /
          n;
    } else {
      static_assert(n >= 2);
      DmPn_of_sin_β[n][m + 1] =
          ((2 * n - 1) * (sin_β * DmPn_of_sin_β[n - 1][m + 1] +
                          (m + 1) * DmPn_of_sin_β[n - 1][m]) -
           (n - 1) * DmPn_of_sin_β[n - 2][m + 1]) /
          n;
    }

#pragma warning(push)
#pragma warning(disable: 4101)
    double cos_β_to_the_m_minus_1;  // Not used if m = 0.
#pragma warning(pop)
    double const 𝔅 = cos_β_to_the_m * DmPn_of_sin_β[n][m];

    double grad_𝔅_polynomials = 0;
    if constexpr (m < n) {
      grad_𝔅_polynomials = cos_β * cos_β_to_the_m * DmPn_of_sin_β[n][m + 1];
    }
    if constexpr (m > 0) {
      cos_β_to_the_m_minus_1 = precomputations.cos_β_to_the_m[m - 1];
      // Remove a singularity when m == 0 and cos_β == 0.
      grad_𝔅_polynomials -=
          m * sin_β * cos_β_to_the_m_minus_1 * DmPn_of_sin_β[n][m];
    }
    Vector<Inverse<Length>, SurfaceFrame> const grad_𝔅 =
        grad_𝔅_polynomials * grad_𝔅_vector;

    double const Cnm = (*cos)[n][m];
    double const Snm = (*sin)[n][m];
    double const 𝔏 = Cnm * cos_mλ + Snm * sin_mλ;

    Vector<Inverse<Length>, SurfaceFrame> 𝔅_grad_𝔏;
    if constexpr (m > 0) {
      // This is not exactly grad_𝔏: we omit the cos_β numerator to remove a
      // singularity.
      Vector<Inverse<Length>, SurfaceFrame> const grad_𝔏 =
          m * (Snm * cos_mλ - Cnm * sin_mλ) * grad_𝔏_vector;
      // Compensate a cos_β to remove a singularity when cos_β == 0.
      𝔅_grad_𝔏 += cos_β_to_the_m_minus_1 * DmPn_of_sin_β[n][m] * grad_𝔏;
    }

    return normalization_factor *
           (grad_ℜ * 𝔅 * 𝔏 + ℜ * grad_𝔅 * 𝔏 + ℜ * 𝔅_grad_𝔏);
  }
}

template<typename Frame>
template<int size, int degree, int... orders>
auto Geopotential<Frame>::
DegreeNAllOrders<size, degree, std::integer_sequence<int, orders...>>::
Acceleration(Displacement<SurfaceFrame> const& r_surface,
             Precomputations<size>& precomputations)
    -> Vector<ReducedAcceleration, SurfaceFrame> {
  if constexpr (degree < 2) {
    return {};
  } else {
    constexpr int n = degree;

    auto const& r² = precomputations.r²;
    auto const& r_norm = precomputations.r_norm;

    auto& ℜ = precomputations.ℜ[n];
    auto& grad_ℜ = precomputations.grad_ℜ;

    // The caller ensures that we process n by increasing values.  Thus, we can
    // safely compute ℜ based on values for lower n's.
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
    grad_ℜ = -(n + 1) * r_surface * ℜ / r²;

    // Force the evaluation by increasing order using an initializer list.
    Accelerations<size> const accelerations = {
        DegreeNOrderM<size, degree, orders>::Acceleration(precomputations)...};

    return (accelerations[orders] + ...);
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
  auto const to_surface_frame = body.ToSurfaceFrame<SurfaceFrame>(t);
  auto const from_surface_frame = to_surface_frame.Inverse();
  Displacement<SurfaceFrame> const r_surface = to_surface_frame(r);

  Precomputations<size> precomputations;

  auto& cos = precomputations.cos;
  auto& sin = precomputations.sin;

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

  auto& DmPn_of_sin_β = precomputations.DmPn_of_sin_β;

  cos = &body.cos();
  sin = &body.sin();

  R3Element<Length> r_surface_coordinates = r_surface.coordinates();
  Length const x = r_surface_coordinates.x;
  Length const y = r_surface_coordinates.y;
  Length const z = r_surface_coordinates.z;

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

  grad_𝔅_vector = UnitVector<SurfaceFrame>({-sin_β * cos_λ,
                                            -sin_β * sin_λ,
                                            cos_β}) / r_norm;
  grad_𝔏_vector = UnitVector<SurfaceFrame>({-sin_λ,
                                            cos_λ,
                                            0}) / r_norm;

  ℜ1 = body.reference_radius() / r²;

  cos_0λ = 1;
  sin_0λ = 0;
  cos_1λ = cos_λ;
  sin_1λ = sin_λ;

  cos_β_to_the_0 = 1;
  cos_β_to_the_1 = cos_β;

  DmPn_of_sin_β[0][0] = 1;
  DmPn_of_sin_β[1][0] = sin_β;
  DmPn_of_sin_β[1][1] = 1;

  // Force the evaluation by increasing degree using an initializer list.
  Accelerations<size> const accelerations = {
      DegreeNAllOrders<size,
                       degrees,
                       std::make_integer_sequence<int, degrees + 1>>::
          Acceleration(r_surface, precomputations)...};

  return from_surface_frame((accelerations[degrees] + ...));
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
  UnitVector<Frame> const& axis = body_->polar_axis();
  return Degree2ZonalAcceleration(axis, r, one_over_r², one_over_r³);
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
    UnitVector<Frame> const& axis,
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

}  // namespace internal_geopotential
}  // namespace physics
}  // namespace principia
