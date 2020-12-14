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


// In the following classes |index| and |indices...| encode the coordinate, the
// degree and the parity (Sin or Cos) of a series.  How this encoding is
// performed determines the order of the series in the basis.

// A helper to build unit polynomials.  |Polynomial| must be a polynomial with
// values in a |dimension|-dimensional Hilbert space.
template<typename Polynomial, int dimension>
struct PolynomialGenerator {
  // Returns a polynomial whose coordinate and degree are encoded in |index|.
  template<int index>
  static Polynomial UnitPolynomial(Instant const& origin);

  // Returns a pair of polynomials (for Sin and Cos) whose parity, coordinate
  // and degree are encoded in |index|.
  template<typename Polynomials, int index>
  static Polynomials UnitPolynomials(Instant const& origin);
};

template<typename Polynomial, int dimension>
template<int index>
Polynomial PolynomialGenerator<Polynomial, dimension>::UnitPolynomial(
    Instant const& origin) {
  // Extract the coordinate and the degree from |index|.  The coordinate varies
  // faster, so terms of the same degree but different coordinates are
  // consecutive when this function is called from a template pack expansion.
  static constexpr int coordinate = index % dimension;
  static constexpr int degree = index / dimension;

  using Coefficients = typename Polynomial::Coefficients;
  using Coefficient = std::tuple_element_t<degree, Coefficients>;
  Coefficients coefficients;
  std::get<degree>(coefficients) =
      CoefficientGenerator<Coefficient, dimension>::template Unit<coordinate>();
  return Polynomial(coefficients, origin);
}

template<typename Polynomial, int dimension>
template<typename Polynomials, int index>
Polynomials PolynomialGenerator<Polynomial, dimension>::UnitPolynomials(
    Instant const& origin) {
  // Extract the parity, coordinate and the degree from |index|.  The coordinate
  // varies faster, so terms of the same parity and degree but different
  // coordinates are consecutive when this function is called from a template
  // pack expansion.  The parity varies next, so Sin and Cos terms alternate for
  // a given degree.  Finally, the degree increases.
  static constexpr int coordinate = index % dimension;
  static constexpr int parity = (index / dimension) % 2;
  static constexpr int degree = index / (2 * dimension);

  // Reencode the degree and coordinate alone for generating the polynomials.
  static constexpr int degree_and_coordinate = coordinate + dimension * degree;

  Polynomial const zero{{}, origin};
  if constexpr (parity == 0) {
    return {/*sin=*/zero,
            /*cos=*/PolynomialGenerator<Polynomial, dimension>::UnitPolynomial<
                degree_and_coordinate>(origin)};
  } else {
    return {/*sin=*/PolynomialGenerator<Polynomial, dimension>::UnitPolynomial<
                degree_and_coordinate>(origin),
            /*cos=*/zero};
  }
}

inline bool PoissonSeriesSubspace::orthogonal(PoissonSeriesSubspace const v,
                                              PoissonSeriesSubspace const w) {
  return v.coordinate_ != w.coordinate_ || v.parity_ != w.parity_;
}

inline PoissonSeriesSubspace::PoissonSeriesSubspace(Coordinate coordinate,
                                                    Parity parity)
    : coordinate_(coordinate), parity_(parity) {}

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
  static std::array<PoissonSeriesSubspace, sizeof...(indices)> Subspaces(
      Instant const& origin);
};

template<typename Series, int degree, int dimension, std::size_t... indices>
std::array<Series, sizeof...(indices)> AperiodicSeriesGenerator<
    Series,
    degree, dimension,
    std::index_sequence<indices...>>::BasisElements(Instant const& origin) {
  return {(Series(
      PolynomialGenerator<typename Series::AperiodicPolynomial,
                          dimension>::template UnitPolynomial<indices>(origin),
      {}))...};
}

template<typename Series, int degree, int dimension, std::size_t... indices>
std::array<PoissonSeriesSubspace, sizeof...(indices)> AperiodicSeriesGenerator<
    Series, degree, dimension,
    std::index_sequence<indices...>>::Subspaces(Instant const& origin) {
  return {PoissonSeriesSubspace{
      static_cast<PoissonSeriesSubspace::Coordinate>(indices % dimension),
      // The degree of the ith polynomial is i / dimension, so its parity is
      // (i / dimension) % 2.
      static_cast<PoissonSeriesSubspace::Parity>((indices / dimension) %
                                                 2)}...};
}


// A helper to generate periodic series, i.e., series of the form  tⁿ sin ω t
// and tⁿ cos ω t along some dimension.  Note how a single index_sequence is
// generated that covers the cross-product of parities, degrees and dimensions.
template<typename Series,
         int degree, int dimension,
         typename = std::make_index_sequence<2 * dimension * (degree + 1)>>
struct PeriodicSeriesGenerator;

template<typename Series, int degree, int dimension, std::size_t... indices>
struct PeriodicSeriesGenerator<Series, degree, dimension,
                               std::index_sequence<indices...>> {
  static std::array<Series, sizeof...(indices)> BasisElements(
      AngularFrequency const& ω,
      Instant const& origin);
  static std::array<PoissonSeriesSubspace, sizeof...(indices)> Subspaces(
      AngularFrequency const& ω,
      Instant const& origin);
};

template<typename Series, int degree, int dimension, std::size_t... indices>
std::array<Series, sizeof...(indices)> PeriodicSeriesGenerator<
    Series,
    degree, dimension,
    std::index_sequence<indices...>>::BasisElements(AngularFrequency const& ω,
                                                    Instant const& origin) {
  typename Series::AperiodicPolynomial const aperiodic_zero{{}, origin};
  return {Series(
      aperiodic_zero,
      {{ω,
        PolynomialGenerator<typename Series::PeriodicPolynomial, dimension>::
            template UnitPolynomials<typename Series::Polynomials, indices>(
                origin)}})...};
}

template<typename Series, int degree, int dimension, std::size_t... indices>
std::array<PoissonSeriesSubspace, sizeof...(indices)> PeriodicSeriesGenerator<
    Series, degree, dimension,
    std::index_sequence<indices...>>::Subspaces(AngularFrequency const& ω,
                                                Instant const& origin) {
  return {
      PoissonSeriesSubspace{
          static_cast<PoissonSeriesSubspace::Coordinate>(indices % dimension),
          // The parity of the trigonometric factors of the ith basis element is
          // ⌊i / dimension⌋ mod 2; the degree of its polynomial factor is
          // ⌊i / (2 * dimension)⌋, so that the overall parity of that basis
          // element is (⌊i / dimension⌋ + ⌊i / (2 * dimension)⌋) mod 2.
          static_cast<PoissonSeriesSubspace::Parity>(
              (indices / dimension + indices / (2 * dimension)) % 2)}...};
}


template<typename Value, int degree,
         template<typename, typename, int> class Evaluator>
auto PoissonSeriesBasisGenerator<Value, degree, Evaluator>::
Basis(Instant const& origin) -> std::array<Series, dimension * (degree + 1)> {
  return AperiodicSeriesGenerator<Series, degree, dimension>::
             BasisElements(origin);
}

template<typename Value, int degree,
         template<typename, typename, int> class Evaluator>
auto PoissonSeriesBasisGenerator<Value, degree, Evaluator>::Subspaces(
    Instant const& origin)
    -> std::array<PoissonSeriesSubspace, dimension*(degree + 1)> {
  return AperiodicSeriesGenerator<Series, degree, dimension>::Subspaces(origin);
}

template<typename Value, int degree,
         template<typename, typename, int> class Evaluator>
auto PoissonSeriesBasisGenerator<Value, degree, Evaluator>::Basis(
    AngularFrequency const& ω,
    Instant const& origin) -> std::array<Series, 2 * dimension * (degree + 1)> {
  return PeriodicSeriesGenerator<Series, degree, dimension>::
             BasisElements(ω, origin);
}

template<typename Value, int degree,
         template<typename, typename, int> class Evaluator>
auto PoissonSeriesBasisGenerator<Value, degree, Evaluator>::Subspaces(
    AngularFrequency const& ω,
    Instant const& origin)
    -> std::array<PoissonSeriesSubspace, 2 * dimension*(degree + 1)> {
  return PeriodicSeriesGenerator<Series, degree, dimension>::Subspaces(ω,
                                                                       origin);
}

inline std::ostream& operator<<(std::ostream& out,
                                PoissonSeriesSubspace const& subspace) {
  return out << "{coordinate: " << static_cast<int>(subspace.coordinate_)
             << ", parity: " << static_cast<int>(subspace.parity_) << "}";
}

}  // namespace internal_poisson_series_basis
}  // namespace numerics
}  // namespace principia
