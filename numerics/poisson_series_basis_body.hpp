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



// A helper to build unit quantities or multivector.  |Coefficient| must be a
// member of a Hilbert space with |dimension| dimensions.
template<typename Coefficient,
         int dimension,
         typename = std::make_index_sequence<dimension>,
         typename = void>
struct CoefficientGenerator;

// Specialization for non-quantities (i.e., Multivector).
template<typename Coefficient, int dimension, std::size_t... ds>
struct CoefficientGenerator<Coefficient,
                            dimension,
                            std::index_sequence<ds...>,
                            std::enable_if_t<!is_quantity_v<Coefficient>>> {
  // Returns a unit multivector in the direction given by |d|, for instance,
  // a multivector with coordinates {1, 0, 0}.
  template<int d>
  static Coefficient Unit();
};

// Specialization for quantities.
template<typename Coefficient, std::size_t... ds>
struct CoefficientGenerator<Coefficient,
                            /*dimension=*/1,
                            std::index_sequence<ds...>,
                            std::enable_if_t<is_quantity_v<Coefficient>>> {
  // Returns a unit quantity.  This function is templated for consistency with
  // the preceeding specialization; |d| must be 0.
  template<int d>
  static Coefficient Unit();
};

template<typename Coefficient, int dimension, std::size_t... ds>
template<int d>
Coefficient
CoefficientGenerator<Coefficient,
                     dimension,
                     std::index_sequence<ds...>,
                     std::enable_if_t<!is_quantity_v<Coefficient>>>::Unit() {
  using Scalar = typename Hilbert<Coefficient>::NormType;
  return Coefficient({(d == ds ? si::Unit<Scalar> : Scalar{})...});
}

template<typename Coefficient, std::size_t... ds>
template<int d>
Coefficient
CoefficientGenerator<Coefficient,
                     /*dimension=*/1,
                     std::index_sequence<ds...>,
                     std::enable_if_t<is_quantity_v<Coefficient>>>::Unit() {
  static_assert(d == 0);
  return si::Unit<Coefficient>;
}


// In the following classes |index| and |indices...| encode the dimension, the
// degree and the parity (Sin or Cos) of a series.  How this encoding is
// performed determines the order of the series in the basis.

// A helper to build unit polynomials.  |Polynomial| must be a polynomial with
// values in a Hilbert space with |dimension| dimensions.
template<typename Polynomial, int dimension>
struct PolynomialGenerator {
  // Returns a polynomial whose dimension and degree are encoded in |index|.
  template<int index>
  static Polynomial UnitPolynomial(Instant const& origin);

  // Returns a pair of polynomials (for Sin and Cos) whose parity, dimension and
  // degree are encoded in |index|.
  template<typename Polynomials, int index>
  static Polynomials UnitPolynomials(Instant const& origin);
};

template<typename Polynomial, int dimension>
template<int index>
Polynomial PolynomialGenerator<Polynomial, dimension>::UnitPolynomial(
    Instant const& origin) {
  // Extract the dimension and the degree from |index|.  The dimension varies
  // faster, so terms of the same degree but different dimensions are
  // consecutive when this function is called from a template pack expansion.
  static constexpr int d = index % dimension;
  static constexpr int n = index / dimension;

  using Coefficients = typename Polynomial::Coefficients;
  using Coefficient = std::tuple_element_t<n, Coefficients>;
  Coefficients coefficients;
  std::get<n>(coefficients) =
      CoefficientGenerator<Coefficient, dimension>::template Unit<d>();
  return Polynomial(coefficients, origin);
}

template<typename Polynomial, int dimension>
template<typename Polynomials, int index>
static Polynomials PolynomialGenerator<Polynomial, dimension>::UnitPolynomials(
    Instant const& origin) {
  // Extract the parity, dimension and the degree from |index|.  The dimension
  // varies faster, so terms of the same parity and degree but different
  // dimensions are consecutive when this function is called from a template
  // pack expansion.  The parity varies next, so Sin and Cos terms alternate for
  // a given degree.  Finally, the degree increases.
  static constexpr int d = index % dimension;
  static constexpr int parity = (index / dimension) % 2;
  static constexpr int n = index / (2 * dimension);

  static Polynomial const zero{{}, origin};
  if constexpr (parity == 0) {
    return {/*sin=*/zero,
            /*cos=*/PolynomialGenerator<Polynomial, dimension>::UnitPolynomial<
                index / 2>(origin)};
  } else {
    return {/*sin=*/PolynomialGenerator<Polynomial, dimension>::UnitPolynomial<
                index / 2>(origin),
            /*cos=*/zero};
  }
}


// A helper to generate aperiodic series, i.e., series of the form tⁿ along some
// dimension.  Note how a single index_sequence is generated that covers the
// cross-product of degrees and dimensions.
template<typename Series,
         int degree, int dimension,
         typename = std::make_index_sequence<dimension * (degree + 1)>>
struct AperiodicSeriesGenerator;

template<typename Series, int degree, int dimension, std::size_t... indices>
struct AperiodicSeriesGenerator<Series,
                                degree, dimension,
                                std::index_sequence<indices...>> {
  static std::array<Series, sizeof...(indices)> BasisElements(
      Instant const& origin);
};

template<typename Series, int degree, int dimension, std::size_t... indices>
std::array<Series, sizeof...(indices)> AperiodicSeriesGenerator<
    Series,
    degree, dimension,
    std::index_sequence<indices...>>::BasisElements(Instant const& origin) {
  return {
      (Series(PolynomialGenerator<typename Series::AperiodicPolynomial,
                                  dimension>::UnitPolynomial<indices>(origin),
              {}))...};
}


// A helper to generate periodic series, i.e., series of the form  tⁿ sin ω t
// and tⁿ cos ω t along some dimension.  Note how a single index_sequence is
// generated that covers the cross-product of parities, degrees and dimensions.
template<typename Series,
         int degree, int dimension,
         typename = std::make_index_sequence<2 * dimension * (degree + 1)>>
struct PeriodicSeriesGenerator;

template<typename Series, int degree, int dimension, std::size_t... indices>
struct PeriodicSeriesGenerator<Series,
                        degree, dimension,
                        std::index_sequence<indices...>> {
  static std::array<Series, sizeof...(indices)> BasisElements(
      AngularFrequency const& ω,
      Instant const& origin);
};

template<typename Series, int degree, int dimension, std::size_t... indices>
std::array<Series, sizeof...(indices)> PeriodicSeriesGenerator<
    Series,
    degree, dimension,
    std::index_sequence<indices...>>::BasisElements(AngularFrequency const& ω,
                                                    Instant const& origin) {
  static typename Series::AperiodicPolynomial const aperiodic_zero{{}, origin};
  return {Series(
      aperiodic_zero,
      {{ω,
        PolynomialGenerator<typename Series::PeriodicPolynomial, dimension>::
            template UnitPolynomials<typename Series::Polynomials, indices>(
                origin)}})...};
}


template<typename Series, int degree>
auto PoissonSeriesBasisGenerator<Series, degree>::Basis(Instant const& origin)
    -> std::array<Series, dimension * (degree + 1)> {
  return AperiodicSeriesGenerator<Series, degree, dimension>::
             BasisElements(origin);
}

template<typename Series, int degree>
auto PoissonSeriesBasisGenerator<Series, degree>::Basis(
    AngularFrequency const& ω,
    Instant const& origin) -> std::array<Series, 2 * dimension * (degree + 1)> {
  return PeriodicSeriesGenerator<Series, degree, dimension>::
             BasisElements(ω, origin);
}

}  // namespace internal_poisson_series_basis
}  // namespace numerics
}  // namespace principia
