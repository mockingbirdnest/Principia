
#pragma once

#include "numerics/frequency_analysis.hpp"

#include <algorithm>
#include <functional>
#include <vector>

#include "base/status.hpp"
#include "base/tags.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/hilbert.hpp"
#include "mathematica/mathematica.hpp"
#include "numerics/poisson_series_basis.hpp"
#include "numerics/root_finders.hpp"
#include "numerics/unbounded_arrays.hpp"
#include "quantities/elementary_functions.hpp"

namespace principia {
namespace numerics {
namespace frequency_analysis {
namespace internal_frequency_analysis {

#define PRINCIPIA_USE_CGS1 0
#define PRINCIPIA_USE_CGS2 0

using base::Error;
using base::Status;
using base::uninitialized;
using geometry::Hilbert;
using geometry::Vector;
using quantities::Inverse;
using quantities::IsFinite;
using quantities::Sqrt;
using quantities::Square;
using quantities::SquareRoot;
using quantities::si::Metre;
using quantities::si::Radian;
using quantities::si::Second;

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

//TODO(phl):comment, parameter names
template<typename BasisSeries,
         int aperiodic_wdegree, int periodic_wdegree,
         template<typename, typename, int> class Evaluator>
Status NormalGramSchmidtStep(
    BasisSeries const& aₘ,
    PoissonSeries<double,
                  aperiodic_wdegree, periodic_wdegree, Evaluator> const& weight,
    Instant const& t_min,
    Instant const& t_max,
    std::vector<PoissonSeriesSubspace> const& basis_subspaces,
    std::vector<BasisSeries> const& q,
    BasisSeries& qₘ,
    UnboundedVector<double>& rₘ) {
  int const m = q.size();
#if PRINCIPIA_USE_CGS1
  //TODO(phl):comment
  // This code follows Björk, Numerics of Gram-Schmidt Orthogonalization,
  // Algorithm 6.1.  It processes one column of q and r at a time.

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
      if (!PoissonSeriesSubspace::orthogonal(basis_subspaces[i],
                                             basis_subspaces[m])) {
        double const sᵖₘ =
            InnerProduct(q[i], previous_q̂ₘ, weight, t_min, t_max);
        q̂ₘ -= sᵖₘ * q[i];
        rₘ[i] += sᵖₘ;
      }
    }
    q̂ₘ_norm = q̂ₘ.Norm(weight, t_min, t_max);

    if (!IsFinite(q̂ₘ_norm)) {
      return Status(Error::OUT_OF_RANGE, u8"Unable to compute q̂ₘ_norm");
    }
    LOG_IF(ERROR, q̂ₘ_norm < α * previous_q̂ₘ_norm)
        << "CGS retrying: " << m << " " << q̂ₘ_norm << " " << previous_q̂ₘ_norm;
  } while (q̂ₘ_norm < α * previous_q̂ₘ_norm);

  // Fill the result.
  qₘ = q̂ₘ / q̂ₘ_norm;
  rₘ[m] = q̂ₘ_norm;
#else
  auto aₘ⁽ᵏ⁾ = aₘ;
  for (int k = 0; k < m; ++k) {
    if (!PoissonSeriesSubspace::orthogonal(basis_subspaces[k],
                                           basis_subspaces[m])) {
      rₘ[k] = InnerProduct(q[k], aₘ⁽ᵏ⁾, weight, t_min, t_max);
      aₘ⁽ᵏ⁾ -= rₘ[k] * q[k];
    }
  }

  auto const rₘₘ = aₘ⁽ᵏ⁾.Norm(weight, t_min, t_max);
  if (!IsFinite(rₘₘ)) {
    return Status(Error::OUT_OF_RANGE, u8"Unable to compute rₘₘ");
  }

  // Fill the result.
  qₘ = aₘ⁽ᵏ⁾ / rₘₘ;
  rₘ[m] = rₘₘ;
#endif
  return Status::OK;
}

//TODO(phl):comment
template<typename Function, typename BasisSeries, typename Norm,
         int aperiodic_wdegree, int periodic_wdegree,
         template<typename, typename, int> class Evaluator>
Status AugmentedGramSchmidtStep(
    Function const& aₘ,
    PoissonSeries<double,
                  aperiodic_wdegree, periodic_wdegree, Evaluator> const& weight,
    Instant const& t_min,
    Instant const& t_max,
    std::vector<BasisSeries> const& q,
    UnboundedVector<Norm>& rₘ) {
  int const m = q.size();
#if PRINCIPIA_USE_CGS2
  //TODO(phl):comment
  // This code follows Björk, Numerics of Gram-Schmidt Orthogonalization,
  // Algorithm 6.1.  It processes one column of q and r at a time.

  static constexpr double α = 0.5;

  Function q̂ₘ = aₘ;
  Norm q̂ₘ_norm = q̂ₘ.Norm(weight, t_min, t_max);

  // Loop on p.
  Function previous_q̂ₘ = q̂ₘ;
  Norm previous_q̂ₘ_norm;
  do {
    previous_q̂ₘ = q̂ₘ;
    previous_q̂ₘ_norm = q̂ₘ_norm;
    for (int i = 0; i < m; ++i) {
      Norm const sᵖₘ = InnerProduct(q[i], previous_q̂ₘ, weight, t_min, t_max);
      q̂ₘ -= sᵖₘ * q[i];
      rₘ[i] += sᵖₘ;
    }
    q̂ₘ_norm = q̂ₘ.Norm(weight, t_min, t_max);

    if (!IsFinite(q̂ₘ_norm)) {
      return Status(Error::OUT_OF_RANGE, u8"Unable to compute q̂ₘ_norm");
    }
  } while (q̂ₘ_norm < α * previous_q̂ₘ_norm);

  // Fill the result.
  rₘ[m] = q̂ₘ_norm;
#else
  auto aₘ⁽ᵏ⁾ = aₘ;
  for (int k = 0; k < m; ++k) {
    rₘ[k] = InnerProduct(q[k], aₘ⁽ᵏ⁾, weight, t_min, t_max);
    aₘ⁽ᵏ⁾ -= rₘ[k] * q[k];
  }

  auto const rₘₘ = aₘ⁽ᵏ⁾.Norm(weight, t_min, t_max);
  if (!IsFinite(rₘₘ)) {
    return Status(Error::OUT_OF_RANGE, u8"Unable to compute rₘₘ");
  }

  // Fill the result.
  rₘ[m] = rₘₘ;
#endif
  return Status::OK;
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
  // step.
  auto const augmented_function = function - F;

  mathematica::Logger logger(TEMP_DIR / "least_square.wl");
  logger.Set(
      "function", function, mathematica::ExpressIn(Metre, Radian, Second));
  logger.Set(
      "tMin", t_min, mathematica::ExpressIn(Metre, Radian, Second));
  logger.Set(
      "tMax", t_max, mathematica::ExpressIn(Metre, Radian, Second));

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

      //LOG_IF(ERROR, m == 0) << qₘ;
      //LOG(ERROR) << "m: " << m << " rₘ[m]: " << rₘ[m]
      //           << " qₘ:" << qₘ.Norm(weight, t_min, t_max);
      for (int i = 0; i <= m; ++i) {
        r[i][m] = rₘ[i];
      }
      q.push_back(qₘ);
      DCHECK_EQ(m + 1, q.size());
    }

    UnboundedVector<Norm> z_ρ(basis_size + 1);
    auto const status = AugmentedGramSchmidtStep(/*aₘ=*/augmented_function,
                                                 weight, t_min, t_max,
                                                 q,
                                                 /*rₘ=*/z_ρ);
    if (!status.ok()) {
      return F;
    }

    logger.Append("q", q, mathematica::ExpressIn(Metre, Radian, Second));
    logger.Append("r", r, mathematica::ExpressIn(Metre, Radian, Second));
    logger.Append("basis", basis, mathematica::ExpressIn(Metre, Radian, Second));
    logger.Append(u8"zρ", z_ρ, mathematica::ExpressIn(Metre, Radian, Second));
    LOG(ERROR) << "ρ: " << z_ρ[z_ρ.size() - 1];
    //for (int i = 0; i < z_ρ.size() - 1; ++i) {
    //  LOG(ERROR) << "z[" << i << "]: " << z_ρ[i];
    //}
    z_ρ.EraseToEnd(z_ρ.size() - 1);
    auto const x = BackSubstitution(r, z_ρ);

    logger.Append("x", x, mathematica::ExpressIn(Metre, Radian, Second));
    F = ResultSeries(result_zero, {{}});
    auto f = augmented_function;
    for (int i = 0; i < x.size(); ++i) {
      auto const x_basis = x[i] * basis[i];
      F += x_basis;
      f -= x_basis;
    }
    logger.Append(
        "approximation", F, mathematica::ExpressIn(Metre, Radian, Second));

    ω = calculator(f);
    if (!ω.has_value()) {
      return F;
    }

    int const ω_basis_size = MakeBasis<aperiodic_degree, periodic_degree>(
        ω, t_min, t_max, basis, basis_subspaces);
    m_begin = basis_size;
    basis_size += ω_basis_size;
    r.Extend(ω_basis_size, uninitialized);
  }
}

}  // namespace internal_frequency_analysis
}  // namespace frequency_analysis
}  // namespace numerics
}  // namespace principia
