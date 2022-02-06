
#pragma once

#include "numerics/frequency_analysis.hpp"

#include <algorithm>
#include <functional>
#include <vector>

#include "absl/status/status.h"
#include "base/tags.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/hilbert.hpp"
#include "numerics/matrix_computations.hpp"
#include "numerics/poisson_series_basis.hpp"
#include "numerics/root_finders.hpp"
#include "numerics/unbounded_arrays.hpp"
#include "quantities/elementary_functions.hpp"

namespace principia {
namespace numerics {
namespace frequency_analysis {
namespace internal_frequency_analysis {

// Note that using CGS and R in some tests produces imprecise results, possibly
// because CGS doesn't yield a good value for R ([Bjö94] is silent on this
// point).  MGS appears a bit faster and more precise than CGS on some tests.
// R has no visible performance effect.
#define PRINCIPIA_USE_CGS 0
#define PRINCIPIA_USE_R 1

using base::uninitialized;
using geometry::Hilbert;
using geometry::Vector;
using quantities::Inverse;
using quantities::IsFinite;
using quantities::Sqrt;
using quantities::Square;
using quantities::SquareRoot;

// Appends basis elements for |ω| to |basis| and |basis_subspaces|.  Returns the
// number of elements that were appended.
template<int aperiodic_degree, int periodic_degree,
         typename BasisSeries>
int MakeBasis(std::optional<AngularFrequency> const& ω,
              Instant const& t_min,
              Instant const& t_max,
              std::vector<BasisSeries>& basis,
              std::vector<PoissonSeriesSubspace>& basis_subspaces) {
  if (ω.value() == AngularFrequency{}) {
    using Generator =
        PoissonSeriesBasisGenerator<BasisSeries, aperiodic_degree>;
    auto const ω_basis = Generator::Basis(t_min, t_max);
    auto const ω_basis_subspaces = Generator::Subspaces();
    std::move(ω_basis.begin(), ω_basis.end(), std::back_inserter(basis));
    std::move(ω_basis_subspaces.begin(),
              ω_basis_subspaces.end(),
              std::back_inserter(basis_subspaces));
    return std::tuple_size_v<decltype(ω_basis)>;
  } else {
    using Generator = PoissonSeriesBasisGenerator<BasisSeries, periodic_degree>;
    auto const ω_basis = Generator::Basis(ω.value(), t_min, t_max);
    auto const ω_basis_subspaces = Generator::Subspaces(ω.value());
    std::move(ω_basis.begin(), ω_basis.end(), std::back_inserter(basis));
    std::move(ω_basis_subspaces.begin(),
              ω_basis_subspaces.end(),
              std::back_inserter(basis_subspaces));
    return std::tuple_size_v<decltype(ω_basis)>;
  }
}

// Given a column |aₘ| of a matrix (or quasimatrix in our case, see [Tre10])
// this function produces the columns |qₘ|, |rₘ| of its QR decomposition.  The
// inner product is defined by |weight|, |t_min| and |t_max|.  |q| is the Q
// quasimatrix constructed so far, and |subspaces| specify the subspaces spanned
// by the |q|s and by |aₘ|.
template<typename BasisSeries,
         int aperiodic_wdegree, int periodic_wdegree,
         template<typename, typename, int> class Evaluator>
absl::Status NormalGramSchmidtStep(
    BasisSeries const& aₘ,
    PoissonSeries<double,
                  aperiodic_wdegree, periodic_wdegree, Evaluator> const& weight,
    Instant const& t_min,
    Instant const& t_max,
    std::vector<PoissonSeriesSubspace> const& subspaces,
    std::vector<BasisSeries> const& q,
    BasisSeries& qₘ,
    UnboundedVector<double>& rₘ) {
  static absl::Status const bad_norm =
    absl::OutOfRangeError("Unable to compute norm");
  int const m = q.size();

#if PRINCIPIA_USE_CGS
  // This code follows [Bjö94], Algorithm 6.1.
  static constexpr double α = 0.5;
  BasisSeries q̂ₘ = aₘ;

  // Formal integration works for a single basis element.
  double q̂ₘ_norm =
      Sqrt((PointwiseInnerProduct(q̂ₘ, q̂ₘ) * weight).Integrate(t_min, t_max) /
           (t_max - t_min));

  // Loop on p.
  BasisSeries previous_q̂ₘ = q̂ₘ;
  double previous_q̂ₘ_norm;
  do {
    previous_q̂ₘ = q̂ₘ;
    previous_q̂ₘ_norm = q̂ₘ_norm;
    for (int i = 0; i < m; ++i) {
      if (!PoissonSeriesSubspace::orthogonal(subspaces[i], subspaces[m])) {
        double const sᵖₘ =
            InnerProduct(q[i], previous_q̂ₘ, weight, t_min, t_max);
        q̂ₘ -= sᵖₘ * q[i];
        rₘ[i] += sᵖₘ;
      }
    }
    q̂ₘ_norm = q̂ₘ.Norm(weight, t_min, t_max);

    if (!IsFinite(q̂ₘ_norm)) {
      return bad_norm;
    }
  } while (q̂ₘ_norm < α * previous_q̂ₘ_norm);

  // Fill the result.
  qₘ = q̂ₘ / q̂ₘ_norm;
  rₘ[m] = q̂ₘ_norm;
#else
  // This code follows [Hig02], Algorithm 19.12.  See also [Bjö94], Algorithm
  // 2.2, for the column version of MGS which is what we are using here.
  auto aₘ⁽ᵏ⁾ = aₘ;
  for (int k = 0; k < m; ++k) {
    if (!PoissonSeriesSubspace::orthogonal(subspaces[k], subspaces[m])) {
      rₘ[k] = InnerProduct(q[k], aₘ⁽ᵏ⁾, weight, t_min, t_max);
      aₘ⁽ᵏ⁾ -= rₘ[k] * q[k];
    }
  }

  auto const rₘₘ = aₘ⁽ᵏ⁾.Norm(weight, t_min, t_max);
  if (!IsFinite(rₘₘ)) {
    return bad_norm;
  }

  // Fill the result.
  qₘ = aₘ⁽ᵏ⁾ / rₘₘ;
  rₘ[m] = rₘₘ;
#endif

  return absl::OkStatus();
}

// This function performs the augmented QR decomposition step described in
// [Hig02] section 20.3.  Note that as an optimization in updates |b|, because
// the computation of |z| for larger and larger R would perform the exact same
// inner products for the range [0, m_begin[.  The range of |q| to process (and
// the range of |z| to update is at indices [m_begin, m_end[.  This function
// doesn't return |qₘ₊₁| because it's not needed for the solution.  It also
// doesn't return |ρ|.
template<typename Function, typename BasisSeries, typename Norm,
         int aperiodic_wdegree, int periodic_wdegree,
         template<typename, typename, int> class Evaluator>
absl::Status AugmentedGramSchmidtStep(
    Function& b,
    PoissonSeries<double,
                  aperiodic_wdegree, periodic_wdegree, Evaluator> const& weight,
    Instant const& t_min,
    Instant const& t_max,
    std::vector<BasisSeries> const& q,
    int const m_begin,
    int const m_end,
    UnboundedVector<Norm>& z) {
  // It would be conceptually possible to use [Bjö94], Algorithm 6.1 here and
  // do reorthonormalization.  Unfortunately, it runs afoul of an issue where
  // the inner product of a piecewise Poisson series with a polynomial doesn't
  // converge (because it depends on a heuristics that uses the maximum
  // frequency).  Instead of trying to make it work, we use MGS.

  // This code follows [Hig02], Algorithm 19.12.  See also [Bjö94], Algorithm
  // 2.2, for the column version of MGS which is what we are using here.
  for (int k = m_begin; k < m_end; ++k) {
    z[k] = InnerProduct(q[k], b, weight, t_min, t_max);
    b -= z[k] * q[k];
  }

  // We do not compute the norm of |b| here (named |ρ| in [Hig02] section 20.3)
  // because it's an additional cost: the client can compute the norm of the
  // residual however they want anyway.

  return absl::OkStatus();
}

template<typename Function,
         int aperiodic_wdegree, int periodic_wdegree,
         template<typename, typename, int> class Evaluator>
AngularFrequency PreciseMode(
    Interval<AngularFrequency> const& fft_mode,
    Function const& function,
    PoissonSeries<double,
                  aperiodic_wdegree, periodic_wdegree,
                  Evaluator> const& weight) {
  auto const weighted_function = weight * function;
  auto const weighted_function_spectrum = weighted_function.FourierTransform();

  auto power =
      [&weighted_function_spectrum](AngularFrequency const& ω) {
        return weighted_function_spectrum(ω).Norm²();
      };

  return Brent(power,
               fft_mode.min,
               fft_mode.max,
               std::greater<>());
}

template<int aperiodic_degree, int periodic_degree,
         typename Function,
         int aperiodic_wdegree, int periodic_wdegree,
         template<typename, typename, int> class Evaluator>
PoissonSeries<std::invoke_result_t<Function, Instant>,
              aperiodic_degree, periodic_degree,
              Evaluator>
Projection(Function const& function,
           AngularFrequency const& ω,
           PoissonSeries<double,
                         aperiodic_wdegree, periodic_wdegree,
                         Evaluator> const& weight,
           Instant const& t_min,
           Instant const& t_max) {
  std::optional<AngularFrequency> optional_ω = ω;

  // A calculator that returns optional_ω once and then stops.
  auto angular_frequency_calculator = [&optional_ω](auto const& residual) {
    auto const result = optional_ω;
    optional_ω = std::nullopt;
    return result;
  };

  return IncrementalProjection<aperiodic_degree, periodic_degree>(
      function,
      angular_frequency_calculator,
      weight,
      t_min, t_max);
}

template<int aperiodic_degree, int periodic_degree,
         typename Function,
         typename AngularFrequencyCalculator,
         int aperiodic_wdegree, int periodic_wdegree,
         template<typename, typename, int> class Evaluator>
PoissonSeries<std::invoke_result_t<Function, Instant>,
              aperiodic_degree, periodic_degree,
              Evaluator>
IncrementalProjection(Function const& function,
                      AngularFrequencyCalculator const& calculator,
                      PoissonSeries<double,
                                    aperiodic_wdegree, periodic_wdegree,
                                    Evaluator> const& weight,
                      Instant const& t_min,
                      Instant const& t_max) {
  using Value = std::invoke_result_t<Function, Instant>;
  using Norm = typename Hilbert<Value>::NormType;
  using Normalized = typename Hilbert<Value>::NormalizedType;
  using BasisSeries = PoissonSeries<Normalized,
                                    aperiodic_degree, periodic_degree,
                                    Evaluator>;
  using ResultSeries = PoissonSeries<Value,
                                     aperiodic_degree, periodic_degree,
                                     Evaluator>;

  Instant const& t0 = weight.origin();
  auto const basis_zero = static_cast<
      typename BasisSeries::AperiodicPolynomial>(
      typename PoissonSeries<Normalized, 0, 0, Evaluator>::AperiodicPolynomial(
          {Normalized{}}, t0));
  auto const result_zero =
      static_cast<typename ResultSeries::AperiodicPolynomial>(
          typename PoissonSeries<Value, 0, 0, Evaluator>::AperiodicPolynomial(
              {Value{}}, t0));

  std::optional<AngularFrequency> ω = calculator(function);
  CHECK(ω.has_value());

  std::vector<BasisSeries> basis;
  // The Poisson series basis[k] belongs to the subspace basis_subspaces[k];
  // this remains true after orthonormalization, i.e., q[k] belongs to the
  // subspace basis_subspaces[k] below.
  std::vector<PoissonSeriesSubspace> basis_subspaces;

  int basis_size = MakeBasis<aperiodic_degree, periodic_degree>(
      ω, t_min, t_max, basis, basis_subspaces);

  // This is logically Q in the QR decomposition of basis.
  std::vector<BasisSeries> q;

  // This is logically R in the QR decomposition of basis.
  UnboundedUpperTriangularMatrix<double> r(basis_size, uninitialized);

  ResultSeries F(result_zero, {{}});

  // The input function with a degree suitable for the augmented Gram-Schmidt
  // step.  Updated by the augmented Gram-Schmidt step.
  auto b = function - F;
  UnboundedVector<Norm> z(basis_size, uninitialized);

  int m_begin = 0;
  for (;;) {
    for (int m = m_begin; m < basis_size; ++m) {
      BasisSeries qₘ(basis_zero, {{}});
      UnboundedVector<double> rₘ(m + 1);

      auto const status = NormalGramSchmidtStep(/*aₘ=*/basis[m],
                                                weight, t_min, t_max,
                                                basis_subspaces, q,
                                                qₘ, rₘ);
      if (!status.ok()) {
        return F;
      }

      // Fill the QR decomposition.
      for (int i = 0; i <= m; ++i) {
        r(i, m) = rₘ[i];
      }
      q.push_back(qₘ);
      DCHECK_EQ(m + 1, q.size());
    }

    auto const status = AugmentedGramSchmidtStep(b,
                                                 weight, t_min, t_max,
                                                 q,
                                                 m_begin, /*m_end=*/basis_size,
                                                 z);
    if (!status.ok()) {
      return F;
    }

    // The conventional way to proceed here ([Hig02], section 20.3, [GV13],
    // section 5.3.5) is to solve R x = z and compute the solution as A x,
    // presumably to get the solution in the canonical basis.  There is no
    // canonical basis for quasimatrices, though, and it's easy to see that the
    // solution can also be expressed as Q z, which appears numerically well-
    // conditioned (note that we don't use R on that path).
#if PRINCIPIA_USE_R
    auto const x = BackSubstitution(r, z);
    F = ResultSeries(result_zero, {{}});
    auto f = function - F;
    for (int i = 0; i < x.size(); ++i) {
      auto const x_basis = x[i] * basis[i];
      F += x_basis;
      f -= x_basis;
    }
#else
    F = ResultSeries(result_zero, {{}});
    auto const f = b;
    for (int i = 0; i < z.size(); ++i) {
      F += z[i] * q[i];
    }
#endif

    ω = calculator(f);
    if (!ω.has_value()) {
      return F;
    }

    int const ω_basis_size = MakeBasis<aperiodic_degree, periodic_degree>(
        ω, t_min, t_max, basis, basis_subspaces);
    m_begin = basis_size;
    basis_size += ω_basis_size;
    r.Extend(ω_basis_size, uninitialized);
    z.Extend(ω_basis_size, uninitialized);
  }
}

#undef PRINCIPIA_USE_CGS
#undef PRINCIPIA_USE_R

}  // namespace internal_frequency_analysis
}  // namespace frequency_analysis
}  // namespace numerics
}  // namespace principia
