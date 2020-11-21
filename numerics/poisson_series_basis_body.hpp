#pragma once

#include "numerics/poisson_series_basis.hpp"

#include <algorithm>

#include "geometry/hilbert.hpp"
#include "quantities/si.hpp"
#include "quantities/traits.hpp"

namespace principia {
namespace numerics {
namespace internal_poisson_series_basis {

using geometry::Hilbert;
using quantities::is_quantity_v;
namespace si = quantities::si;

// A helper struct for generating the Poisson series tⁿ sin ω t and tⁿ cos ω t.
// d is either 0 (for a 1-dimensional value type) or 0, 1, 2 (for a
// 3-dimensional value type).
template<typename Series, int n, int d>
struct SeriesGenerator {
  // The series tⁿ.
  static Series Aperiodic(Instant const& origin);
  // The series tⁿ sin ω t.
  static Series Sin(AngularFrequency const& ω, Instant const& origin);
  // The series tⁿ cos ω t.
  static Series Cos(AngularFrequency const& ω, Instant const& origin);

 private:
  // The polynomial tⁿ.
  template<typename Polynomial>
  static Polynomial Unit(Instant const& origin);
};

template<typename Series, int n, int d>
Series SeriesGenerator<Series, n, d>::Aperiodic(Instant const& origin) {
  return Series(Unit<typename Series::AperiodicPolynomial>(origin), {});
}

template<typename Series, int n, int d>
Series SeriesGenerator<Series, n, d>::Sin(AngularFrequency const& ω,
                                          Instant const& origin) {
  typename Series::AperiodicPolynomial const aperiodic_zero{{}, origin};
  typename Series::PeriodicPolynomial const periodic_zero{{}, origin};
  return Series(aperiodic_zero,
                {{ω,
                  {/*sin=*/Unit<typename Series::PeriodicPolynomial>(origin),
                   /*cos=*/periodic_zero}}});
}

template<typename Series, int n, int d>
Series SeriesGenerator<Series, n, d>::Cos(AngularFrequency const& ω,
                                          Instant const& origin) {
  typename Series::AperiodicPolynomial const aperiodic_zero{{}, origin};
  typename Series::PeriodicPolynomial const periodic_zero{{}, origin};
  return Series(
      aperiodic_zero,
      {{ω,
        {/*sin=*/periodic_zero,
         /*cos=*/Unit<typename Series::PeriodicPolynomial>(origin)}}});
}

template<typename Series, int n, int d>
template<typename Polynomial>
Polynomial SeriesGenerator<Series, n, d>::Unit(Instant const& origin) {
  typename Polynomial::Coefficients coefficients;
  using Coefficient =
      std::tuple_element_t<n, typename Polynomial::Coefficients>;
  Coefficient& coefficient = std::get<n>(coefficients);
  if constexpr (is_quantity_v<Coefficient>) {
    coefficient = si::Unit<Coefficient>;
  } else {
    using Scalar = typename Hilbert<Coefficient>::NormType;
    if constexpr (d == 0) {
      coefficient = Coefficient({si::Unit<Scalar>, Scalar{}, Scalar{}});
    } else if constexpr (d == 1) {
      coefficient = Coefficient({Scalar{}, si::Unit<Scalar>, Scalar{}});
    } else if constexpr (d == 2) {
      coefficient = Coefficient({Scalar{}, Scalar{}, si::Unit<Scalar>});
    }
  }
  return Polynomial(coefficients, origin);
}



template<typename Coefficient,
         int dimensions,
         typename = std::make_index_sequence<dimensions>,
         typename = void>
struct CoefficientGenerator;

template<typename Coefficient, int dimensions, std::size_t... ds>
struct CoefficientGenerator<Coefficient,
                            dimensions,
                            std::index_sequence<ds...>,
                            std::enable_if_t<!is_quantity_v<Coefficient>>> {
  template<int d>
  static Coefficient Unit() {
    using Scalar = typename Hilbert<Coefficient>::NormType;
    return Coefficient({(d == ds ? si::Unit<Scalar> : Scalar{})...});
  }
};

template<typename Coefficient, std::size_t... ds>
struct CoefficientGenerator<Coefficient,
                            1,
                            std::index_sequence<ds...>,
                            std::enable_if_t<is_quantity_v<Coefficient>>> {
  template<int d>
  static Coefficient Unit() {
    static_assert(d == 0);
    return si::Unit<Coefficient>;
  }
};




template<typename Polynomial, int dimensions>
struct PolynomialGenerator {
  template<int index>
  static Polynomial Unit(Instant const& origin) {
    static constexpr int d = index % dimensions;
    static constexpr int n = index / dimensions;
    using Coefficients = typename Polynomial::Coefficients;
    using Coefficient = std::tuple_element_t<n, Coefficients>;
    Coefficients coefficients;
    std::get<n>(coefficients) =
        CoefficientGenerator<Coefficient, dimensions>::template Unit<d>();
    return Polynomial(coefficients, origin);
  }
};


template<typename Series,
         int degree,
         int dimensions,
         typename = std::make_index_sequence<dimensions * (degree + 1)>>
struct SeriesGenerator2;

template<typename Series, int degree, int dimensions, std::size_t... indices>
struct SeriesGenerator2<Series,
                        degree,
                        dimensions,
                        std::index_sequence<indices...>> {
  static std::array<Series, sizeof...(indices)> Make(Instant const& origin) {
    return {(Series(PolynomialGenerator<typename Series::AperiodicPolynomial,
                                        dimensions>::Unit<indices>(origin),
                    {}))...};
  }
};



template<typename Series,
         int degree,
         int dimensions,
         typename = std::make_index_sequence<2 * dimensions * (degree + 1)>>
struct SeriesGenerator3;

template<typename Series, int degree, int dimensions, std::size_t... indices>
struct SeriesGenerator3<Series,
                        degree,
                        dimensions,
                        std::index_sequence<indices...>> {
  template<int index>
  static typename Series::Polynomials Unit(Instant const& origin) {
    static constexpr int d = index % dimensions;
    static constexpr int parity = (index / dimensions) % 2;
    static constexpr int n = index / (2 * dimensions);
    static typename Series::PeriodicPolynomial const zero{{}, origin};
    if constexpr (parity == 0) {
      return {/*sin=*/zero,
              /*cos=*/PolynomialGenerator<typename Series::PeriodicPolynomial,
                                          dimensions>::Unit<index / 2>(origin)};
    } else {
      return {/*sin=*/PolynomialGenerator<typename Series::PeriodicPolynomial,
                                          dimensions>::Unit<index / 2>(origin),
              /*cos=*/zero};
    }
  }

  static std::array<Series, sizeof...(indices)> Make(AngularFrequency const& ω,
                                                     Instant const& origin) {
    static
        typename Series::AperiodicPolynomial const aperiodic_zero{{}, origin};
    return {Series(aperiodic_zero, {{ω, Unit<indices>(origin)}})...};
  }
};



template<typename Series, int degree>
auto PoissonSeriesBasisGenerator<Series, degree>::Basis(Instant const& origin)
    -> std::array<Series, dimension * (degree + 1)> {
  return SeriesGenerator2<Series, degree, dimension>::Make(origin);
}

template<typename Series, int degree>
auto PoissonSeriesBasisGenerator<Series, degree>::Basis(
    AngularFrequency const& ω,
    Instant const& origin) -> std::array<Series, 2 * dimension * (degree + 1)> {
  return SeriesGenerator3<Series, degree, dimension>::Make(ω, origin);
}

}  // namespace internal_poisson_series_basis
}  // namespace numerics
}  // namespace principia
