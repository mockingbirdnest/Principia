#pragma once

#include "physics/geopotential.hpp"

#include <algorithm>
#include <cmath>
#include <queue>
#include <vector>

#include "base/tags.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/r3_element.hpp"
#include "numerics/fixed_arrays.hpp"
#include "numerics/legendre_normalization_factor.mathematica.h"
#include "numerics/max_abs_normalized_associated_legendre_function.mathematica.h"
#include "numerics/polynomial_evaluators.hpp"
#include "quantities/elementary_functions.hpp"
#include "quantities/quantities.hpp"

namespace principia {
namespace physics {
namespace internal_geopotential {

using numerics::FixedLowerTriangularMatrix;
using numerics::FixedVector;
using numerics::HornerEvaluator;
using numerics::LegendreNormalizationFactor;
using numerics::MaxAbsNormalizedAssociatedLegendreFunction;
using namespace principia::base::_tags;
using namespace principia::geometry::_grassmann;
using namespace principia::geometry::_r3_element;
using namespace principia::quantities::_elementary_functions;
using namespace principia::quantities::_named_quantities;
using namespace principia::quantities::_quantities;

// The notation in this file follows documentation/Geopotential.pdf.

template<typename Frame>
struct Geopotential<Frame>::Precomputations {
  // Allocate the maximum size to cover all possible degrees.  Making |size| a
  // template parameter of this class would be possible, but it would greatly
  // increase the number of instances of DegreeNOrderM and friends.
  static constexpr int size = OblateBody<Frame>::max_geopotential_degree + 1;

  // These quantities are independent from n and m.
  typename OblateBody<Frame>::GeopotentialCoefficients const* cos;
  typename OblateBody<Frame>::GeopotentialCoefficients const* sin;

  Length r_norm;
  Square<Length> r²;
  Vector<double, Frame> r_normalized;  // Only used for the acceleration.

  double sin_β;
  double cos_β;

  // Only used for the acceleration.
  Vector<double, Frame> grad_𝔅_vector;
  Vector<double, Frame> grad_𝔏_vector;

  // These quantities depend on n but are independent from m.
  FixedVector<Exponentiation<Length, -2>, size> ℜ_over_r{
      uninitialized};  // 0 unused.

  // These quantities depend on m but are independent from n.
  FixedVector<double, size> cos_mλ{uninitialized};  // 0 unused.
  FixedVector<double, size> sin_mλ{uninitialized};  // 0 unused.
  FixedVector<double, size> cos_β_to_the_m{uninitialized};

  // These quantities depend on both n and m.  Note that the zeros for m > n are
  // not stored.
  FixedLowerTriangularMatrix<double, size> DmPn_of_sin_β{uninitialized};
};

template<typename Frame>
template<int degree, int order>
class Geopotential<Frame>::DegreeNOrderM {
 public:
  static auto Acceleration(
      Inverse<Square<Length>> const& σℜ_over_r,
      Vector<Inverse<Square<Length>>, Frame> const& grad_σℜ,
      Precomputations& precomputations) -> Vector<ReducedAcceleration, Frame>;

  static auto Potential(Inverse<Square<Length>> const& σℜ_over_r,
                        Precomputations& precomputations) -> ReducedPotential;

 private:
  static void UpdatePrecomputations(Precomputations& precomputations);
};

template<typename Frame>
template<int degree, int... orders>
class Geopotential<Frame>::
DegreeNAllOrders<degree, std::integer_sequence<int, orders...>> {
 public:
  static auto Acceleration(Geopotential<Frame> const& geopotential,
                           Precomputations& precomputations)
      -> Vector<ReducedAcceleration, Frame>;

  static auto Potential(Geopotential<Frame> const& geopotential,
                        Precomputations& precomputations) -> ReducedPotential;

 private:
  static void UpdatePrecomputations(Precomputations& precomputations);
};

template<typename Frame>
template<int... degrees>
class Geopotential<Frame>::AllDegrees<std::integer_sequence<int, degrees...>> {
 public:
  static auto Acceleration(Geopotential<Frame> const& geopotential,
                           Instant const& t,
                           Displacement<Frame> const& r,
                           Length const& r_norm,
                           Square<Length> const& r²,
                           Exponentiation<Length, -3> const& one_over_r³)
      -> Vector<ReducedAcceleration, Frame>;

  static auto Potential(Geopotential<Frame> const& geopotential,
                        Instant const& t,
                        Displacement<Frame> const& r,
                        Length const& r_norm,
                        Square<Length> const& r²,
                        Exponentiation<Length, -3> const& one_over_r³)
      -> ReducedPotential;

 private:
  static void InitializePrecomputations(
      Geopotential<Frame> const& geopotential,
      Instant const& t,
      Displacement<Frame> const& r,
      Length const& r_norm,
      Square<Length> const& r²,
      Exponentiation<Length, -3> const& one_over_r³,
      Precomputations& precomputations);
};

template<typename Frame>
template<int degree, int order>
auto Geopotential<Frame>::DegreeNOrderM<degree, order>::Acceleration(
    Inverse<Square<Length>> const& σℜ_over_r,
    Vector<Inverse<Square<Length>>, Frame> const& grad_σℜ,
    Precomputations& precomputations)
    -> Vector<ReducedAcceleration, Frame> {
  UpdatePrecomputations(precomputations);

  if constexpr (degree == 2 && order == 1) {
    return {};
  } else {
    constexpr int n = degree;
    constexpr int m = order;
    static_assert(0 <= m && m <= n);

    double const cos_β = precomputations.cos_β;
    double const sin_β = precomputations.sin_β;

    auto const& grad_𝔅_vector = precomputations.grad_𝔅_vector;
    auto const& grad_𝔏_vector = precomputations.grad_𝔏_vector;

    auto const& cos_mλ = precomputations.cos_mλ[m];
    auto const& sin_mλ = precomputations.sin_mλ[m];

    auto const& cos_β_to_the_m = precomputations.cos_β_to_the_m[m];

    auto const& DmPn_of_sin_β = precomputations.DmPn_of_sin_β;
    auto const& cos = *precomputations.cos;
    auto const& sin = *precomputations.sin;

    constexpr double normalization_factor = LegendreNormalizationFactor(n, m);

#pragma warning(push)
#pragma warning(disable: 4101)
    double cos_β_to_the_m_minus_1;  // Not used if m = 0.
#pragma warning(pop)
    double const 𝔅 = cos_β_to_the_m * DmPn_of_sin_β(n, m);

    double grad_𝔅_polynomials = 0;
    if constexpr (m < n) {
      grad_𝔅_polynomials = cos_β * cos_β_to_the_m * DmPn_of_sin_β(n, m + 1);
    }
    if constexpr (m > 0) {
      cos_β_to_the_m_minus_1 = precomputations.cos_β_to_the_m[m - 1];
      // Remove a singularity when m == 0 and cos_β == 0.
      grad_𝔅_polynomials -=
          m * sin_β * cos_β_to_the_m_minus_1 * DmPn_of_sin_β(n, m);
    }

    double const Cnm = cos(n, m);
    double const Snm = sin(n, m);
    double 𝔏;
    if constexpr (m == 0) {
      𝔏 = Cnm;
    } else {
      𝔏 = Cnm * cos_mλ + Snm * sin_mλ;
    }

    Vector<ReducedAcceleration, Frame> const 𝔅𝔏_grad_ℜ = (𝔅 * 𝔏) * grad_σℜ;
    Vector<ReducedAcceleration, Frame> const ℜ𝔏_grad_𝔅 =
        (σℜ_over_r * 𝔏 * grad_𝔅_polynomials) * grad_𝔅_vector;
    Vector<ReducedAcceleration, Frame> grad_ℜ𝔅𝔏 = 𝔅𝔏_grad_ℜ + ℜ𝔏_grad_𝔅;
    if constexpr (m > 0) {
      // Compensate a cos_β to remove a singularity when cos_β == 0.
      Vector<ReducedAcceleration, Frame> const ℜ𝔅_grad_𝔏 =
          (σℜ_over_r *
           cos_β_to_the_m_minus_1 * DmPn_of_sin_β(n, m) *  // 𝔅/cos_β
           m * (Snm * cos_mλ - Cnm * sin_mλ)) * grad_𝔏_vector;  // grad_𝔏*cos_β
      grad_ℜ𝔅𝔏 += ℜ𝔅_grad_𝔏;
    }

    return normalization_factor * grad_ℜ𝔅𝔏;
  }
}

template<typename Frame>
template<int degree, int order>
auto Geopotential<Frame>::DegreeNOrderM<degree, order>::Potential(
    Inverse<Square<Length>> const& σℜ_over_r,
    Precomputations& precomputations) -> ReducedPotential {
  UpdatePrecomputations(precomputations);

  if constexpr (degree == 2 && order == 1) {
    return ReducedPotential{};
  } else {
    constexpr int n = degree;
    constexpr int m = order;
    static_assert(0 <= m && m <= n);

    auto const& r_norm = precomputations.r_norm;

    auto const& cos_mλ = precomputations.cos_mλ[m];
    auto const& sin_mλ = precomputations.sin_mλ[m];

    auto const& cos_β_to_the_m = precomputations.cos_β_to_the_m[m];

    auto const& DmPn_of_sin_β = precomputations.DmPn_of_sin_β;
    auto const& cos = *precomputations.cos;
    auto const& sin = *precomputations.sin;

    constexpr double normalization_factor = LegendreNormalizationFactor(n, m);

    Inverse<Length> const σℜ = r_norm * σℜ_over_r;
    double const 𝔅 = cos_β_to_the_m * DmPn_of_sin_β(n, m);

    double const Cnm = cos(n, m);
    double const Snm = sin(n, m);
    double 𝔏;
    if constexpr (m == 0) {
      𝔏 = Cnm;
    } else {
      𝔏 = Cnm * cos_mλ + Snm * sin_mλ;
    }

    return -normalization_factor * σℜ * 𝔅 * 𝔏;
  }
}

template<typename Frame>
template<int degree, int order>
void Geopotential<Frame>::DegreeNOrderM<degree, order>::UpdatePrecomputations(
    Precomputations& precomputations) {
  if constexpr (degree == 2 && order == 1) {
    // Let's not forget the Legendre derivative that we would compute if we did
    // not short-circuit.
    precomputations.DmPn_of_sin_β(2, 2) = 3;
  } else {
    constexpr int n = degree;
    constexpr int m = order;
    static_assert(0 <= m && m <= n);

    double const sin_β = precomputations.sin_β;

    auto& cos_mλ = precomputations.cos_mλ[m];
    auto& sin_mλ = precomputations.sin_mλ[m];

    auto& cos_β_to_the_m = precomputations.cos_β_to_the_m[m];
    auto& DmPn_of_sin_β = precomputations.DmPn_of_sin_β;

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
        cos_mλ = (cos_hλ + sin_hλ) * (cos_hλ - sin_hλ);
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
      DmPn_of_sin_β(n, 0) = ((2 * n - 1) * sin_β * DmPn_of_sin_β(n - 1, 0) -
                             (n - 1) * DmPn_of_sin_β(n - 2, 0)) /
                            n;
    }

    // Recurrence relationship between the associated Legendre polynomials.
    // Account for the fact that DmPn_of_sin_β is identically zero if m > n.
    if constexpr (m == n) {
      // Do not store the zero.
    } else if constexpr (m == n - 1) {  // NOLINT(readability/braces)
      static_assert(n >= 1);
      DmPn_of_sin_β(n, m + 1) =
          ((2 * n - 1) * (m + 1) * DmPn_of_sin_β(n - 1, m)) / n;
    } else if constexpr (m == n - 2) {  // NOLINT(readability/braces)
      static_assert(n >= 1);
      DmPn_of_sin_β(n, m + 1) =
          ((2 * n - 1) * (sin_β * DmPn_of_sin_β(n - 1, m + 1) +
                          (m + 1) * DmPn_of_sin_β(n - 1, m))) /
          n;
    } else {
      static_assert(n >= 2);
      DmPn_of_sin_β(n, m + 1) =
          ((2 * n - 1) * (sin_β * DmPn_of_sin_β(n - 1, m + 1) +
                          (m + 1) * DmPn_of_sin_β(n - 1, m)) -
           (n - 1) * DmPn_of_sin_β(n - 2, m + 1)) /
          n;
    }
  }
}

template<typename Frame>
template<int degree, int... orders>
auto Geopotential<Frame>::
DegreeNAllOrders<degree, std::integer_sequence<int, orders...>>::
Acceleration(Geopotential<Frame> const& geopotential,
             Precomputations& precomputations)
    -> Vector<ReducedAcceleration, Frame> {
  if constexpr (degree < 2) {
    return {};
  } else {
    constexpr int n = degree;
    constexpr int size = sizeof...(orders);

    UpdatePrecomputations(precomputations);

    auto const& r_norm = precomputations.r_norm;
    auto const& r² = precomputations.r²;
    auto const& r_normalized = precomputations.r_normalized;

    auto const& ℜ_over_r = precomputations.ℜ_over_r[n];
    auto const ℜʹ = -(n + 1) * ℜ_over_r;
    // Note that ∇ℜ = ℜʹ * r_normalized.

    Inverse<Square<Length>> σℜ_over_r;
    Vector<Inverse<Square<Length>>, Frame> grad_σℜ;
    if constexpr (n == 2 && size > 1) {
      geopotential.degree_damping_[2].ComputeDampedRadialQuantities(
          r_norm,
          r²,
          r_normalized,
          ℜ_over_r,
          ℜʹ,
          σℜ_over_r,
          grad_σℜ);
      // If we are above the outer threshold, we should not have been called
      // (σ = 0).
      DCHECK_LT(r_norm, geopotential.degree_damping_[2].outer_threshold());
      Vector<ReducedAcceleration, Frame> const j2_acceleration =
          DegreeNOrderM<2, 0>::Acceleration(
              σℜ_over_r, grad_σℜ, precomputations);
      geopotential.sectoral_damping_.ComputeDampedRadialQuantities(
          r_norm,
          r²,
          r_normalized,
          ℜ_over_r,
          ℜʹ,
          σℜ_over_r,
          grad_σℜ);
      // If we are above the outer threshold, we should have been called with
      // (orders...) = (0).
      DCHECK_LT(r_norm, geopotential.sectoral_damping_.outer_threshold());
      // Perform the precomputations for order 1 (but the result is known to be
      // 0, so don't bother adding it).
      DegreeNOrderM<2, 1>::Acceleration(
          σℜ_over_r, grad_σℜ, precomputations);
      Vector<ReducedAcceleration, Frame> const c22_s22_acceleration =
          DegreeNOrderM<2, 2>::Acceleration(
              σℜ_over_r, grad_σℜ, precomputations);
      return j2_acceleration + c22_s22_acceleration;
    } else {
      geopotential.degree_damping_[n].ComputeDampedRadialQuantities(
          r_norm,
          r²,
          r_normalized,
          ℜ_over_r,
          ℜʹ,
          σℜ_over_r,
          grad_σℜ);
      // If we are above the outer threshold, we should not have been called
      // (σ = 0).
      DCHECK_LT(r_norm, geopotential.degree_damping_[n].outer_threshold());

      // Force the evaluation by increasing order using an initializer list.
      ReducedAccelerations<size> const accelerations = {
          DegreeNOrderM<degree, orders>::Acceleration(
              σℜ_over_r, grad_σℜ, precomputations)...};

      return (accelerations[orders] + ...);
    }
  }
}

template<typename Frame>
template<int degree, int... orders>
auto Geopotential<Frame>::
DegreeNAllOrders<degree, std::integer_sequence<int, orders...>>::
Potential(Geopotential<Frame> const& geopotential,
          Precomputations& precomputations) -> ReducedPotential {
  if constexpr (degree < 2) {
    return ReducedPotential{};
  } else {
    constexpr int n = degree;
    constexpr int size = sizeof...(orders);

    UpdatePrecomputations(precomputations);

    auto const& r_norm = precomputations.r_norm;
    auto const& r² = precomputations.r²;
    auto const& ℜ_over_r = precomputations.ℜ_over_r[n];

    Inverse<Square<Length>> σℜ_over_r;
    if constexpr (n == 2 && size > 1) {
      geopotential.degree_damping_[2].ComputeDampedRadialQuantities(r_norm,
                                                                    r²,
                                                                    ℜ_over_r,
                                                                    σℜ_over_r);
      // If we are above the outer threshold, we should not have been called
      // (σ = 0).
      DCHECK_LT(r_norm, geopotential.degree_damping_[2].outer_threshold());
      ReducedPotential const j2_potential =
          DegreeNOrderM<2, 0>::Potential(σℜ_over_r, precomputations);
      geopotential.sectoral_damping_.ComputeDampedRadialQuantities(r_norm,
                                                                   r²,
                                                                   ℜ_over_r,
                                                                   σℜ_over_r);
      // If we are above the outer threshold, we should have been called with
      // (orders...) = (0).
      DCHECK_LT(r_norm, geopotential.sectoral_damping_.outer_threshold());
      // Perform the precomputations for order 1 (but the result is known to be
      // 0, so don't bother adding it).
      DegreeNOrderM<2, 1>::Potential(σℜ_over_r, precomputations);
      ReducedPotential const c22_s22_potential =
          DegreeNOrderM<2, 2>::Potential(σℜ_over_r, precomputations);
      return j2_potential + c22_s22_potential;
    } else {
      geopotential.degree_damping_[n].ComputeDampedRadialQuantities(r_norm,
                                                                    r²,
                                                                    ℜ_over_r,
                                                                    σℜ_over_r);
      // If we are above the outer threshold, we should not have been called
      // (σ = 0).
      DCHECK_LT(r_norm, geopotential.degree_damping_[n].outer_threshold());

      // Force the evaluation by increasing order using an initializer list.
      ReducedPotentials<size> const potentials = {
          DegreeNOrderM<degree, orders>::Potential(σℜ_over_r,
                                                   precomputations)...};

      return (potentials[orders] + ...);
    }
  }
}

template<typename Frame>
template<int degree, int... orders>
void Geopotential<Frame>::
DegreeNAllOrders<degree, std::integer_sequence<int, orders...>>::
UpdatePrecomputations(Precomputations& precomputations) {
  constexpr int n = degree;

  auto const& r² = precomputations.r²;
  auto& ℜ_over_r = precomputations.ℜ_over_r[n];

  // The caller ensures that we process n by increasing values.  Thus, we can
  // safely compute ℜ based on values for lower n's.
  if constexpr (n % 2 == 0) {
    int const h = n / 2;
    auto const& ℜh_over_r = precomputations.ℜ_over_r[h];
    ℜ_over_r = ℜh_over_r * ℜh_over_r * r²;
  } else {
    int const h1 = n / 2;
    int const h2 = n - h1;
    auto const& ℜh1_over_r = precomputations.ℜ_over_r[h1];
    auto const& ℜh2_over_r = precomputations.ℜ_over_r[h2];
    ℜ_over_r = ℜh1_over_r * ℜh2_over_r * r²;
  }
}

template<typename Frame>
template<int... degrees>
auto Geopotential<Frame>::AllDegrees<std::integer_sequence<int, degrees...>>::
Acceleration(Geopotential<Frame> const& geopotential,
             Instant const& t,
             Displacement<Frame> const& r,
             Length const& r_norm,
             Square<Length> const& r²,
             Exponentiation<Length, -3> const& one_over_r³)
    -> Vector<ReducedAcceleration, Frame> {
  constexpr int size = sizeof...(degrees);
  OblateBody<Frame> const& body = *geopotential.body_;
  const bool is_zonal =
      body.is_zonal() ||
      r_norm > geopotential.sectoral_damping_.outer_threshold();

  Precomputations precomputations;
  InitializePrecomputations(
      geopotential, t, r, r_norm, r², one_over_r³, precomputations);

  // Force the evaluation by increasing degree using an initializer list.  In
  // the zonal case, no point in going beyond order 0.
  ReducedAccelerations<size> accelerations;
  if (is_zonal) {
    accelerations = {
        DegreeNAllOrders<degrees, std::make_integer_sequence<int, 1>>::
            Acceleration(geopotential, precomputations)...};
  } else {
    accelerations = {
        DegreeNAllOrders<degrees,
                         std::make_integer_sequence<int, degrees + 1>>::
            Acceleration(geopotential, precomputations)...};
  }

  return (accelerations[degrees] + ...);
}

template<typename Frame>
template<int... degrees>
auto Geopotential<Frame>::AllDegrees<std::integer_sequence<int, degrees...>>::
Potential(Geopotential<Frame> const& geopotential,
          Instant const& t,
          Displacement<Frame> const& r,
          Length const& r_norm,
          Square<Length> const& r²,
          Exponentiation<Length, -3> const& one_over_r³)
    -> ReducedPotential {
  constexpr int size = sizeof...(degrees);
  OblateBody<Frame> const& body = *geopotential.body_;
  const bool is_zonal =
      body.is_zonal() ||
      r_norm > geopotential.sectoral_damping_.outer_threshold();

  Precomputations precomputations;
  InitializePrecomputations(
      geopotential, t, r, r_norm, r², one_over_r³, precomputations);

  // Force the evaluation by increasing degree using an initializer list.  In
  // the zonal case, no point in going beyond order 0.
  ReducedPotentials<size> potentials;
  if (is_zonal) {
    potentials = {
        DegreeNAllOrders<degrees, std::make_integer_sequence<int, 1>>::
            Potential(geopotential, precomputations)...};
  } else {
    potentials = {
        DegreeNAllOrders<degrees,
                         std::make_integer_sequence<int, degrees + 1>>::
            Potential(geopotential, precomputations)...};
  }

  return (potentials[degrees] + ...);
}

template<typename Frame>
template<int... degrees>
void Geopotential<Frame>::AllDegrees<std::integer_sequence<int, degrees...>>::
InitializePrecomputations(Geopotential<Frame> const& geopotential,
                          Instant const& t,
                          Displacement<Frame> const& r,
                          Length const& r_norm,
                          Square<Length> const& r²,
                          Exponentiation<Length, -3> const& one_over_r³,
                          Precomputations& precomputations) {
  OblateBody<Frame> const& body = *geopotential.body_;
  const bool is_zonal =
      body.is_zonal() ||
      r_norm > geopotential.sectoral_damping_.outer_threshold();

  precomputations.r_norm = r_norm;
  precomputations.r² = r²;

  auto& cos = precomputations.cos;
  auto& sin = precomputations.sin;

  auto& cos_β = precomputations.cos_β;
  auto& sin_β = precomputations.sin_β;

  auto& grad_𝔅_vector = precomputations.grad_𝔅_vector;
  auto& grad_𝔏_vector = precomputations.grad_𝔏_vector;

  auto& ℜ1_over_r = precomputations.ℜ_over_r[1];

  auto& cos_1λ = precomputations.cos_mλ[1];
  auto& sin_1λ = precomputations.sin_mλ[1];

  auto& cos_β_to_the_0 = precomputations.cos_β_to_the_m[0];
  auto& cos_β_to_the_1 = precomputations.cos_β_to_the_m[1];

  auto& DmPn_of_sin_β = precomputations.DmPn_of_sin_β;

  // In the zonal case the rotation of the body is of no importance, so any pair
  // of equatorial vectors will do.
  UnitVector x̂;
  UnitVector ŷ;
  UnitVector const ẑ = body.polar_axis();
  if (is_zonal) {
    x̂ = body.equatorial();
    ŷ = body.biequatorial();
  } else {
    auto const from_surface_frame =
      body.template FromSurfaceFrame<SurfaceFrame>(t);
    x̂ = from_surface_frame(x_);
    ŷ = from_surface_frame(y_);
  }

  Length const x = InnerProduct(r, x̂);
  Length const y = InnerProduct(r, ŷ);
  Length const z = InnerProduct(r, ẑ);

  Square<Length> const x²_plus_y² = x * x + y * y;
  Length const r_equatorial = Sqrt(x²_plus_y²);

  // TODO(phl): This is probably incorrect for celestials that don't have
  // longitudes counted to the East.
  double cos_λ = 1;
  double sin_λ = 0;
  if (r_equatorial > Length{}) {
    Inverse<Length> const one_over_r_equatorial = 1 / r_equatorial;
    cos_λ = x * one_over_r_equatorial;
    sin_λ = y * one_over_r_equatorial;
  }

  cos = &body.cos();
  sin = &body.sin();

  Inverse<Length> const one_over_r_norm = 1 / r_norm;
  precomputations.r_normalized = r * one_over_r_norm;

  cos_β = r_equatorial * one_over_r_norm;
  sin_β = z * one_over_r_norm;

  grad_𝔅_vector = (-sin_β * cos_λ) * x̂ - (sin_β * sin_λ) * ŷ + cos_β * ẑ;
  grad_𝔏_vector = cos_λ * ŷ - sin_λ * x̂;

  ℜ1_over_r = body.reference_radius() * one_over_r³;

  cos_1λ = cos_λ;
  sin_1λ = sin_λ;

  cos_β_to_the_0 = 1;
  cos_β_to_the_1 = cos_β;

  DmPn_of_sin_β(0, 0) = 1;
  DmPn_of_sin_β(1, 0) = sin_β;
  DmPn_of_sin_β(1, 1) = 1;
}

template<typename Frame>
Geopotential<Frame>::Geopotential(not_null<OblateBody<Frame> const*> body,
                                  double const tolerance)
    : body_(body) {
  CHECK_GE(tolerance, 0);
  double const& ε = tolerance;

  // Thresholds for individual harmonics, with lexicographic (threshold, order,
  // degree) comparison.
  // Note that the order of the fields is (degree, order) as usual; comparison
  // is order-first in the priority queue.
  struct Threshold {
    Length r;
    int n;
    int m;
  };
  // If |after(left, right)|, |left| is popped after |right| in the
  // |priority_queue|.
  auto const after = [](Threshold const& left, Threshold const& right) -> bool {
    return left.r < right.r ||
           (left.r == right.r &&
            (left.m < right.m || (left.m == right.m && left.n < right.n)));
  };
  std::priority_queue<Threshold, std::vector<Threshold>, decltype(after)>
      harmonic_thresholds(after);
  for (int n = 2; n <= body_->geopotential_degree(); ++n) {
    for (int m = 0; m <= n; ++m) {
      double const max_abs_Pnm =
          MaxAbsNormalizedAssociatedLegendreFunction(n, m);
      double const Cnm = body->cos()(n, m);
      double const Snm = body->sin()(n, m);
      // TODO(egg): write a rootn.
      Length const r = Cnm == 0 && Snm == 0
                           ? Length{}
                           : body->reference_radius() *
                                 std::pow((max_abs_Pnm * (n + 1) *
                                           Sqrt(Pow<2>(Cnm) + Pow<2>(Snm))) /
                                              ε,
                                          1.0 / n);
      harmonic_thresholds.push({r, n, m});
    }
  }

  harmonic_thresholds.push({Infinity<Length>, 0, 0});
  harmonic_thresholds.push({Infinity<Length>, 1, 0});

  while (!harmonic_thresholds.empty()) {
    auto const& threshold = harmonic_thresholds.top();
    if (threshold.n == 2 && threshold.m == 2) {
      if (degree_damping_.size() > 3) {
        // Enforce the monotonicity relation for sectoral damping.
        sectoral_damping_ =
            HarmonicDamping(degree_damping_[3].inner_threshold());
      } else {
        sectoral_damping_ = HarmonicDamping(threshold.r);
      }
    }
    // Make the thresholds monotonic, using the degree n threshold for all
    // degrees k < n that would otherwise have a lower threshold.
    while (threshold.n >= degree_damping_.size()) {
      degree_damping_.emplace_back(threshold.r);
    }
    harmonic_thresholds.pop();
  }
}

#define PRINCIPIA_CASE_SPHERICAL_HARMONICS_ACCELERATION(d)                     \
  case (d):                                                                    \
    return AllDegrees<std::make_integer_sequence<int, (d) + 1>>::Acceleration( \
        *this, t, r, r_norm, r², one_over_r³)

template<typename Frame>
Vector<Quotient<Acceleration, GravitationalParameter>, Frame>
Geopotential<Frame>::GeneralSphericalHarmonicsAcceleration(
    Instant const& t,
    Displacement<Frame> const& r,
    Length const& r_norm,
    Square<Length> const& r²,
    Exponentiation<Length, -3> const& one_over_r³) const {
  if (r_norm != r_norm) {
    // Short-circuit NaN, to avoid having to deal with an unordered
    // |r_norm| when finding the partition point below.
    return NaN<ReducedAcceleration> * Vector<double, Frame>{};
  }
  // We have |max_degree > 0|.
  int const max_degree = LimitingDegree(r_norm) - 1;
  switch (max_degree) {
    PRINCIPIA_CASE_SPHERICAL_HARMONICS_ACCELERATION(2);
    PRINCIPIA_CASE_SPHERICAL_HARMONICS_ACCELERATION(3);
    PRINCIPIA_CASE_SPHERICAL_HARMONICS_ACCELERATION(4);
    PRINCIPIA_CASE_SPHERICAL_HARMONICS_ACCELERATION(5);
    PRINCIPIA_CASE_SPHERICAL_HARMONICS_ACCELERATION(6);
    PRINCIPIA_CASE_SPHERICAL_HARMONICS_ACCELERATION(7);
    PRINCIPIA_CASE_SPHERICAL_HARMONICS_ACCELERATION(8);
    PRINCIPIA_CASE_SPHERICAL_HARMONICS_ACCELERATION(9);
    PRINCIPIA_CASE_SPHERICAL_HARMONICS_ACCELERATION(10);
    PRINCIPIA_CASE_SPHERICAL_HARMONICS_ACCELERATION(11);
    PRINCIPIA_CASE_SPHERICAL_HARMONICS_ACCELERATION(12);
    PRINCIPIA_CASE_SPHERICAL_HARMONICS_ACCELERATION(13);
    PRINCIPIA_CASE_SPHERICAL_HARMONICS_ACCELERATION(14);
    PRINCIPIA_CASE_SPHERICAL_HARMONICS_ACCELERATION(15);
    PRINCIPIA_CASE_SPHERICAL_HARMONICS_ACCELERATION(16);
    PRINCIPIA_CASE_SPHERICAL_HARMONICS_ACCELERATION(17);
    PRINCIPIA_CASE_SPHERICAL_HARMONICS_ACCELERATION(18);
    PRINCIPIA_CASE_SPHERICAL_HARMONICS_ACCELERATION(19);
    PRINCIPIA_CASE_SPHERICAL_HARMONICS_ACCELERATION(20);
    PRINCIPIA_CASE_SPHERICAL_HARMONICS_ACCELERATION(21);
    PRINCIPIA_CASE_SPHERICAL_HARMONICS_ACCELERATION(22);
    PRINCIPIA_CASE_SPHERICAL_HARMONICS_ACCELERATION(23);
    PRINCIPIA_CASE_SPHERICAL_HARMONICS_ACCELERATION(24);
    PRINCIPIA_CASE_SPHERICAL_HARMONICS_ACCELERATION(25);
    PRINCIPIA_CASE_SPHERICAL_HARMONICS_ACCELERATION(26);
    PRINCIPIA_CASE_SPHERICAL_HARMONICS_ACCELERATION(27);
    PRINCIPIA_CASE_SPHERICAL_HARMONICS_ACCELERATION(28);
    PRINCIPIA_CASE_SPHERICAL_HARMONICS_ACCELERATION(29);
    PRINCIPIA_CASE_SPHERICAL_HARMONICS_ACCELERATION(30);
#if PRINCIPIA_GEOPOTENTIAL_MAX_DEGREE_50
    PRINCIPIA_CASE_SPHERICAL_HARMONICS_ACCELERATION(31);
    PRINCIPIA_CASE_SPHERICAL_HARMONICS_ACCELERATION(32);
    PRINCIPIA_CASE_SPHERICAL_HARMONICS_ACCELERATION(33);
    PRINCIPIA_CASE_SPHERICAL_HARMONICS_ACCELERATION(34);
    PRINCIPIA_CASE_SPHERICAL_HARMONICS_ACCELERATION(35);
    PRINCIPIA_CASE_SPHERICAL_HARMONICS_ACCELERATION(36);
    PRINCIPIA_CASE_SPHERICAL_HARMONICS_ACCELERATION(37);
    PRINCIPIA_CASE_SPHERICAL_HARMONICS_ACCELERATION(38);
    PRINCIPIA_CASE_SPHERICAL_HARMONICS_ACCELERATION(39);
    PRINCIPIA_CASE_SPHERICAL_HARMONICS_ACCELERATION(40);
    PRINCIPIA_CASE_SPHERICAL_HARMONICS_ACCELERATION(41);
    PRINCIPIA_CASE_SPHERICAL_HARMONICS_ACCELERATION(42);
    PRINCIPIA_CASE_SPHERICAL_HARMONICS_ACCELERATION(43);
    PRINCIPIA_CASE_SPHERICAL_HARMONICS_ACCELERATION(44);
    PRINCIPIA_CASE_SPHERICAL_HARMONICS_ACCELERATION(45);
    PRINCIPIA_CASE_SPHERICAL_HARMONICS_ACCELERATION(46);
    PRINCIPIA_CASE_SPHERICAL_HARMONICS_ACCELERATION(47);
    PRINCIPIA_CASE_SPHERICAL_HARMONICS_ACCELERATION(48);
    PRINCIPIA_CASE_SPHERICAL_HARMONICS_ACCELERATION(49);
    PRINCIPIA_CASE_SPHERICAL_HARMONICS_ACCELERATION(50);
#endif
    case 1:
      return Vector<ReducedAcceleration, Frame>{};
    default:
      LOG(FATAL) << "Unexpected degree " << max_degree << " " << body_->name();
      base::noreturn();
  }
}

#undef PRINCIPIA_CASE_SPHERICAL_HARMONICS_ACCELERATION

#define PRINCIPIA_CASE_SPHERICAL_HARMONICS_POTENTIAL(d)                     \
  case (d):                                                                 \
    return AllDegrees<std::make_integer_sequence<int, (d) + 1>>::Potential( \
        *this, t, r, r_norm, r², one_over_r³)

template<typename Frame>
Quotient<SpecificEnergy, GravitationalParameter>
Geopotential<Frame>::GeneralSphericalHarmonicsPotential(
    Instant const& t,
    Displacement<Frame> const& r,
    Length const& r_norm,
    Square<Length> const& r²,
    Exponentiation<Length, -3> const& one_over_r³) const {
  if (r_norm != r_norm) {
    // Short-circuit NaN, to avoid having to deal with an unordered
    // |r_norm| when finding the partition point below.
    return NaN<ReducedPotential>;
  }
  // We have |max_degree > 0|.
  int const max_degree = LimitingDegree(r_norm) - 1;
  switch (max_degree) {
    PRINCIPIA_CASE_SPHERICAL_HARMONICS_POTENTIAL(2);
    PRINCIPIA_CASE_SPHERICAL_HARMONICS_POTENTIAL(3);
    PRINCIPIA_CASE_SPHERICAL_HARMONICS_POTENTIAL(4);
    PRINCIPIA_CASE_SPHERICAL_HARMONICS_POTENTIAL(5);
    PRINCIPIA_CASE_SPHERICAL_HARMONICS_POTENTIAL(6);
    PRINCIPIA_CASE_SPHERICAL_HARMONICS_POTENTIAL(7);
    PRINCIPIA_CASE_SPHERICAL_HARMONICS_POTENTIAL(8);
    PRINCIPIA_CASE_SPHERICAL_HARMONICS_POTENTIAL(9);
    PRINCIPIA_CASE_SPHERICAL_HARMONICS_POTENTIAL(10);
    PRINCIPIA_CASE_SPHERICAL_HARMONICS_POTENTIAL(11);
    PRINCIPIA_CASE_SPHERICAL_HARMONICS_POTENTIAL(12);
    PRINCIPIA_CASE_SPHERICAL_HARMONICS_POTENTIAL(13);
    PRINCIPIA_CASE_SPHERICAL_HARMONICS_POTENTIAL(14);
    PRINCIPIA_CASE_SPHERICAL_HARMONICS_POTENTIAL(15);
    PRINCIPIA_CASE_SPHERICAL_HARMONICS_POTENTIAL(16);
    PRINCIPIA_CASE_SPHERICAL_HARMONICS_POTENTIAL(17);
    PRINCIPIA_CASE_SPHERICAL_HARMONICS_POTENTIAL(18);
    PRINCIPIA_CASE_SPHERICAL_HARMONICS_POTENTIAL(19);
    PRINCIPIA_CASE_SPHERICAL_HARMONICS_POTENTIAL(20);
    PRINCIPIA_CASE_SPHERICAL_HARMONICS_POTENTIAL(21);
    PRINCIPIA_CASE_SPHERICAL_HARMONICS_POTENTIAL(22);
    PRINCIPIA_CASE_SPHERICAL_HARMONICS_POTENTIAL(23);
    PRINCIPIA_CASE_SPHERICAL_HARMONICS_POTENTIAL(24);
    PRINCIPIA_CASE_SPHERICAL_HARMONICS_POTENTIAL(25);
    PRINCIPIA_CASE_SPHERICAL_HARMONICS_POTENTIAL(26);
    PRINCIPIA_CASE_SPHERICAL_HARMONICS_POTENTIAL(27);
    PRINCIPIA_CASE_SPHERICAL_HARMONICS_POTENTIAL(28);
    PRINCIPIA_CASE_SPHERICAL_HARMONICS_POTENTIAL(29);
    PRINCIPIA_CASE_SPHERICAL_HARMONICS_POTENTIAL(30);
#if PRINCIPIA_GEOPOTENTIAL_MAX_DEGREE_50
    PRINCIPIA_CASE_SPHERICAL_HARMONICS_POTENTIAL(31);
    PRINCIPIA_CASE_SPHERICAL_HARMONICS_POTENTIAL(32);
    PRINCIPIA_CASE_SPHERICAL_HARMONICS_POTENTIAL(33);
    PRINCIPIA_CASE_SPHERICAL_HARMONICS_POTENTIAL(34);
    PRINCIPIA_CASE_SPHERICAL_HARMONICS_POTENTIAL(35);
    PRINCIPIA_CASE_SPHERICAL_HARMONICS_POTENTIAL(36);
    PRINCIPIA_CASE_SPHERICAL_HARMONICS_POTENTIAL(37);
    PRINCIPIA_CASE_SPHERICAL_HARMONICS_POTENTIAL(38);
    PRINCIPIA_CASE_SPHERICAL_HARMONICS_POTENTIAL(39);
    PRINCIPIA_CASE_SPHERICAL_HARMONICS_POTENTIAL(40);
    PRINCIPIA_CASE_SPHERICAL_HARMONICS_POTENTIAL(41);
    PRINCIPIA_CASE_SPHERICAL_HARMONICS_POTENTIAL(42);
    PRINCIPIA_CASE_SPHERICAL_HARMONICS_POTENTIAL(43);
    PRINCIPIA_CASE_SPHERICAL_HARMONICS_POTENTIAL(44);
    PRINCIPIA_CASE_SPHERICAL_HARMONICS_POTENTIAL(45);
    PRINCIPIA_CASE_SPHERICAL_HARMONICS_POTENTIAL(46);
    PRINCIPIA_CASE_SPHERICAL_HARMONICS_POTENTIAL(47);
    PRINCIPIA_CASE_SPHERICAL_HARMONICS_POTENTIAL(48);
    PRINCIPIA_CASE_SPHERICAL_HARMONICS_POTENTIAL(49);
    PRINCIPIA_CASE_SPHERICAL_HARMONICS_POTENTIAL(50);
#endif
    case 1:
      return ReducedPotential{};
    default:
      LOG(FATAL) << "Unexpected degree " << max_degree << " " << body_->name();
      base::noreturn();
  }
}

#undef PRINCIPIA_CASE_SPHERICAL_HARMONICS_POTENTIAL

template<typename Frame>
std::vector<HarmonicDamping> const& Geopotential<Frame>::degree_damping()
    const {
  return degree_damping_;
}

template<typename Frame>
HarmonicDamping const& Geopotential<Frame>::sectoral_damping() const {
  return sectoral_damping_;
}

template<typename Frame>
int Geopotential<Frame>::LimitingDegree(Length const& r_norm) const {
  return std::partition_point(
             degree_damping_.begin(),
             degree_damping_.end(),
             [r_norm](HarmonicDamping const& degree_damping) -> bool {
               return r_norm < degree_damping.outer_threshold();
             }) -
         degree_damping_.begin();
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
