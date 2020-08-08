
#pragma once

#include "numerics/frequency_analysis.hpp"

#include <functional>

#include "numerics/root_finders.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace numerics {
namespace frequency_analysis {
namespace internal_frequency_analysis {

using quantities::Square;
namespace si = quantities::si;

template<typename Function,
         typename RValue, int rdegree_, int wdegree_,
         template<typename, typename, int> class Evaluator>
AngularFrequency PreciseMode(
    Interval<AngularFrequency> const& fft_mode,
    Function const& function,
    PoissonSeries<double, wdegree_, Evaluator> const& weight,
    DotProduct<Function, RValue, rdegree_, wdegree_, Evaluator> const& dot) {
  using DotResult =
      Primitive<Product<std::invoke_result_t<Function, Instant>, RValue>, Time>;
  using Degree0 = PoissonSeries<double, 0, Evaluator>;

  auto amplitude = [&dot, &function, &weight](AngularFrequency const& ω) {
    Instant const& t0 = weight.origin();
    Degree0 const sin(typename Degree0::Polynomial({0}, t0),
                      {{ω,
                        {/*sin=*/typename Degree0::Polynomial({1}, t0),
                         /*cos=*/typename Degree0::Polynomial({0}, t0)}}});
    Degree0 const cos(typename Degree0::Polynomial({0}, t0),
                      {{ω,
                        {/*sin=*/typename Degree0::Polynomial({0}, t0),
                         /*cos=*/typename Degree0::Polynomial({1}, t0)}}});
    auto const sin_amplitude = dot(function, sin, weight);
    auto const cos_amplitude = dot(function, cos, weight);
    return sin_amplitude * sin_amplitude + cos_amplitude * cos_amplitude;
  };

  return GoldenSectionSearch(amplitude,
                             fft_mode.min,
                             fft_mode.max,
                             std::greater<Square<DotResult>>());
}

template<typename Polynomial, int index>
Polynomial One(Instant const& origin) {
  typename Polynomial::Coefficients coefficients;
  std::get<index>(coefficients) =
      si::Unit<std::tuple_element_t<index, typename Polynomial::Coefficients>>;
  return Polynomial(coefficients, origin);
}

template<typename Series, int index>
struct SeriesGenerator {
  static Series Sin(AngularFrequency const& ω, Instant const& origin) {
    typename Series::Polynomial::Coefficients const zeros;
    typename Series::Polynomial const zero{zeros, origin};
    return Series(zero,
                  {{ω,
                    {/*sin=*/One<typename Series::Polynomial, index>(origin),
                     /*cos=*/zero}}});
  }
  static Series Cos(AngularFrequency const& ω, Instant const& origin) {
    typename Series::Polynomial::Coefficients const zeros;
    typename Series::Polynomial const zero{zeros, origin};
    return Series(
        zero,
        {{ω,
          {/*sin=*/zero,
           /*cos=*/One<typename Series::Polynomial, index>(origin)}}});
  }
};

template<typename Series,
         typename = std::make_index_sequence<Series::degree>>
struct BasisGenerator;

template<typename Series, std::size_t... indices>
struct BasisGenerator<Series, std::index_sequence<indices...>> {
  static std::array<Series, 2 * Series::degree + 1> Basis(
      AngularFrequency const& ω,
      Instant const& origin) {
    // This has the elements {1, t Sin(ωt), t² Sin(ωt), ... t Cos(ωt)...}
    // which is not the order we want (we want lower-degree polynomials first).
    std::array<Series, 2 * Series::degree + 1> all_series = {
        Series(One<typename Series::Polynomial, 0>(origin), {{}}),
        SeriesGenerator<Series, indices>::Sin(ω, origin)...,
        SeriesGenerator<Series, indices>::Cos(ω, origin)...};
    if (Series::degree >= 2) {
      // Reads the series at index i, compute the index j where it must be
      // stored, squirrel out the series at index j, store the series that was
      // at index i, repeat.
      //TODO(phl):Do we believe that there is a single cycle????
      int i = 2;
      auto si = all_series[i];
      do {
        int const j =
            i <= Series::degree ? 2 * i - 1 : 2 * (i - Series::degree);
        auto const sj = all_series[j];
        all_series[j] = si;
        i = j;
        si = sj;
      } while (i != 2);
    }
    return all_series;
  }
};

template<typename Function,
         typename RValue, int rdegree_, int wdegree_,
         template<typename, typename, int> class Evaluator>
PoissonSeries<std::invoke_result_t<Function, Instant>, rdegree_, Evaluator>
Projection(
    AngularFrequency const& ω,
    Function const& function,
    PoissonSeries<double, wdegree_, Evaluator> const& weight,
    DotProduct<Function, RValue, rdegree_, wdegree_, Evaluator> const& dot) {
  using Value = std::invoke_result_t<Function, Instant>;
  using Series = PoissonSeries<Value, degree_, Evaluator>;
  Instant const& t0 = weight.origin();

  std::vector<Series> basis;
  typename Series::Polynomial::Coefficients const zeros;
  basis.emplace_back(typename Series::Polynomial({1}, t0), {{}});
  typename Series::Polynomial const zero({0}, t0);
  for (int d = 1; d <= degree_; ++d) {
    basis.emplace_back(zero,
                      {{ω,
                        {/*sin=*/typename Series::Polynomial({1}, t0),
                         /*cos=*/typename Series::Polynomial({0}, t0)}}});
    basis.emplace_back(zero,
                      {{ω,
                        {/*sin=*/typename Series::Polynomial({1}, t0),
                         /*cos=*/typename Series::Polynomial({0}, t0)}}});
  }
}

}  // namespace internal_frequency_analysis
}  // namespace frequency_analysis
}  // namespace numerics
}  // namespace principia
