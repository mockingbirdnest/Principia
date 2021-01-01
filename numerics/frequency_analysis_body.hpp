
#pragma once

#include "numerics/frequency_analysis.hpp"

#include <algorithm>
#include <functional>
#include <vector>

#include "base/tags.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/hilbert.hpp"
#include "numerics/poisson_series_basis.hpp"
#include "numerics/root_finders.hpp"
#include "numerics/unbounded_arrays.hpp"
#include "quantities/elementary_functions.hpp"

namespace principia {
namespace numerics {
namespace frequency_analysis {
namespace internal_frequency_analysis {

using base::uninitialized;
using geometry::Hilbert;
using geometry::Vector;
using quantities::Inverse;
using quantities::IsFinite;
using quantities::Sqrt;
using quantities::Square;
using quantities::SquareRoot;

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

  // This code follows [Kud07], section 2.  Our indices start at 0, unlike those
  // of Кудрявцев which start at 1.

  Instant const& t0 = weight.origin();

  std::optional<AngularFrequency> ω = calculator(function);
  CHECK(ω.has_value());

  std::vector<BasisSeries> basis;
  // The Poisson series basis[k] belongs to the subspace basis_subspaces[k];
  // this remains true after orthonormalization, i.e., q[k] belongs to the
  // subspace basis_subspaces[k] below.
  std::vector<PoissonSeriesSubspace> basis_subspaces;

  int basis_size;
  // TODO(phl): This is replicated below.
  if (ω.value() == AngularFrequency{}) {
    using Generator = PoissonSeriesBasisGenerator<BasisSeries,
                                                  aperiodic_degree>;
    auto const ω_basis = Generator::Basis(t_min, t_max);
    auto const ω_basis_subspaces = Generator::Subspaces();
    basis_size = std::tuple_size_v<decltype(ω_basis)>;
    std::move(ω_basis.begin(), ω_basis.end(), std::back_inserter(basis));
    std::move(ω_basis_subspaces.begin(),
              ω_basis_subspaces.end(),
              std::back_inserter(basis_subspaces));
  } else {
    using Generator = PoissonSeriesBasisGenerator<BasisSeries,
                                                  periodic_degree>;
    auto const ω_basis = Generator::Basis(ω.value(), t_min, t_max);
    auto const ω_basis_subspaces = Generator::Subspaces(ω.value());
    basis_size = std::tuple_size_v<decltype(ω_basis)>;
    std::move(ω_basis.begin(), ω_basis.end(), std::back_inserter(basis));
    std::move(ω_basis_subspaces.begin(),
              ω_basis_subspaces.end(),
              std::back_inserter(basis_subspaces));
  }

  // This is logically Q in the QR decomposition of basis.
  std::vector<BasisSeries> q;

  auto const& a₀ = basis[0];
  auto const r₀₀ = a₀.Norm(weight, t_min, t_max);
  CHECK(IsFinite(r₀₀)) << a₀;
  q.push_back(a₀ / r₀₀);

  auto const A₀ = InnerProduct(function, q[0], weight, t_min, t_max);

  ResultSeries F = A₀ * q[0];
  auto f = function - F;

  int m_begin = 1;
  for (;;) {
    for (int m = m_begin; m < basis_size; ++m) {
      auto const& aₘ = basis[m];

      // This code follows Björk, Numerics of Gram-Schmidt Orthogonalization,
      // Algorithm 6.1.
      static constexpr double α = 0.5;
      BasisSeries q̂ₘ = aₘ;

      // Formal integration works for a single basis element.
      double q̂ₘ_norm = Sqrt(
          (PointwiseInnerProduct(q̂ₘ, q̂ₘ) * weight).Integrate(t_min, t_max) /
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
            auto const sᵖₘ =
                InnerProduct(q[i], previous_q̂ₘ, weight, t_min, t_max);
            q̂ₘ -= sᵖₘ * q[i];
          }
        }
        q̂ₘ_norm = q̂ₘ.Norm(weight, t_min, t_max);

        if (!IsFinite(q̂ₘ_norm)) {
          return F;
        }
      } while (q̂ₘ_norm < α * previous_q̂ₘ_norm);

      q.push_back(q̂ₘ / q̂ₘ_norm);
      DCHECK_EQ(m + 1, q.size());

      Norm const Aₘ = InnerProduct(f, q[m], weight, t_min, t_max);

      auto const Aₘqₘ = Aₘ * q[m];
      f -= Aₘqₘ;
      F += Aₘqₘ;
    }

    ω = calculator(f);
    if (!ω.has_value()) {
      return F;
    }

    int ω_basis_size;
    if (ω.value() == AngularFrequency{}) {
      using Generator = PoissonSeriesBasisGenerator<BasisSeries,
                                                    aperiodic_degree>;
      auto const ω_basis = Generator::Basis(t_min, t_max);
      auto const ω_basis_subspaces = Generator::Subspaces();
      ω_basis_size = std::tuple_size_v<decltype(ω_basis)>;
      std::move(ω_basis.begin(), ω_basis.end(), std::back_inserter(basis));
      std::move(ω_basis_subspaces.begin(),
                ω_basis_subspaces.end(),
                std::back_inserter(basis_subspaces));
    } else {
      using Generator = PoissonSeriesBasisGenerator<BasisSeries,
                                                    periodic_degree>;
      auto const ω_basis = Generator::Basis(ω.value(), t_min, t_max);
      auto const ω_basis_subspaces = Generator::Subspaces(ω.value());
      ω_basis_size = std::tuple_size_v<decltype(ω_basis)>;
      std::move(ω_basis.begin(), ω_basis.end(), std::back_inserter(basis));
      std::move(ω_basis_subspaces.begin(),
                ω_basis_subspaces.end(),
                std::back_inserter(basis_subspaces));
    }
    m_begin = basis_size;
    basis_size += ω_basis_size;
  }
}

}  // namespace internal_frequency_analysis
}  // namespace frequency_analysis
}  // namespace numerics
}  // namespace principia
