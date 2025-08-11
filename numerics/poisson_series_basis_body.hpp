#pragma once

#include "numerics/poisson_series_basis.hpp"

#include <algorithm>

#include "geometry/barycentre_calculator.hpp"
#include "numerics/elementary_functions.hpp"

namespace principia {
namespace numerics {
namespace _poisson_series_basis {
namespace internal {

using namespace principia::geometry::_barycentre_calculator;
using namespace principia::numerics::_elementary_functions;


// A helper to build unit quantities or multivector.  `Coefficient` must be a
// member of a Hilbert space with `dimension` dimensions, and must be free from
// Quantity.
template<typename Coefficient,
         int dimension,
         typename = std::make_index_sequence<dimension>,
         typename = void>
struct CoefficientGenerator;

// Specialization for non-double (i.e., Multivector<double, ...>).
template<typename Coefficient, int dimension, std::size_t... ds>
struct CoefficientGenerator<Coefficient,
                            dimension,
                            std::index_sequence<ds...>> {
  // Returns a unit multivector in the direction given by `d`, for instance,
  // a multivector with coordinates {1, 0, 0}.
  template<int d>
  static Coefficient Unit();
};

// Specialization for double.
template<std::size_t... ds>
struct CoefficientGenerator<double,
                            /*dimension=*/1,
                            std::index_sequence<ds...>> {
  // Returns a unit quantity.  This function is templated for consistency with
  // the preceeding specialization; `d` must be 0.
  template<int d>
  static double Unit();
};

template<typename Coefficient, int dimension, std::size_t... ds>
template<int d>
Coefficient CoefficientGenerator<Coefficient,
                                 dimension,
                                 std::index_sequence<ds...>>::Unit() {
  return Coefficient({(d == ds ? 1 : 0)...});
}

template<std::size_t... ds>
template<int d>
double CoefficientGenerator<double,
                            /*dimension=*/1,
                            std::index_sequence<ds...>>::Unit() {
  static_assert(d == 0);
  return 1;
}


// In the following classes `index` and `indices...` encode the coordinate, the
// degree and the parity (Sin or Cos) of a series.  How this encoding is
// performed determines the order of the series in the basis.

// A helper to build unit polynomials.  `Polynomial` must be a polynomial with
// values in a `dimension`-dimensional Hilbert space.
template<typename Polynomial, int dimension>
struct PolynomialGenerator {
  // Returns a polynomial whose coordinate and degree are encoded in `index`.
  template<int index>
  static Polynomial UnitPolynomial(Instant const& t_min,
                                   Instant const& t_mid,
                                   Instant const& t_max);

  // Returns a pair of polynomials (for Sin and Cos) whose parity, coordinate
  // and degree are encoded in `index`.
  template<typename Polynomials, int index>
  static Polynomials UnitPolynomials(Instant const& t_min,
                                     Instant const& t_mid,
                                     Instant const& t_max);
};

template<typename Polynomial, int dimension>
template<int index>
Polynomial PolynomialGenerator<Polynomial, dimension>::UnitPolynomial(
    Instant const& t_min,
    Instant const& t_mid,
    Instant const& t_max) {
  // Extract the coordinate and the degree from `index`.  The coordinate varies
  // faster, so terms of the same degree but different coordinates are
  // consecutive when this function is called from a template pack expansion.
  static constexpr int coordinate = index % dimension;
  static constexpr int degree = index / dimension;

  using Coefficients = typename Polynomial::Coefficients;
  using Coefficient = std::tuple_element_t<degree, Coefficients>;
  using NormalizedCoefficient = typename Hilbert<Coefficient>::NormalizedType;
  Coefficients coefficients;
  std::get<degree>(coefficients) =
      CoefficientGenerator<NormalizedCoefficient,
                           dimension>::template Unit<coordinate>() /
      Pow<degree>(0.5 * (t_max - t_min));
  return Polynomial(coefficients, t_mid);
}

template<typename Polynomial, int dimension>
template<typename Polynomials, int index>
Polynomials PolynomialGenerator<Polynomial, dimension>::UnitPolynomials(
    Instant const& t_min,
    Instant const& t_mid,
    Instant const& t_max) {
  // Extract the parity, coordinate and the degree from `index`.  The coordinate
  // varies faster, so terms of the same parity and degree but different
  // coordinates are consecutive when this function is called from a template
  // pack expansion.  The parity varies next, so Sin and Cos terms alternate for
  // a given degree.  Finally, the degree increases.
  static constexpr int coordinate = index % dimension;
  static constexpr int parity = (index / dimension) % 2;
  static constexpr int degree = index / (2 * dimension);

  // Reencode the degree and coordinate alone for generating the polynomials.
  static constexpr int degree_and_coordinate = coordinate + dimension * degree;

  Polynomial const zero{{}, t_mid};
  if constexpr (parity == 0) {
    return {.sin = zero,
            .cos = PolynomialGenerator<Polynomial, dimension>::UnitPolynomial<
                degree_and_coordinate>(t_min, t_mid, t_max)};
  } else {
    return {.sin = PolynomialGenerator<Polynomial, dimension>::UnitPolynomial<
                degree_and_coordinate>(t_min, t_mid, t_max),
            .cos = zero};
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
      Instant const& t_min,
      Instant const& t_max);
  static std::array<PoissonSeriesSubspace, sizeof...(indices)> Subspaces();
};

template<typename Series, int degree, int dimension, std::size_t... indices>
std::array<Series, sizeof...(indices)> AperiodicSeriesGenerator<
    Series,
    degree, dimension,
    std::index_sequence<indices...>>::BasisElements(Instant const& t_min,
                                                    Instant const& t_max) {
  Instant const t_mid = Barycentre({t_min, t_max});
  return {(Series(
      PolynomialGenerator<typename Series::AperiodicPolynomial,
                          dimension>::template UnitPolynomial<indices>(t_min,
                                                                       t_mid,
                                                                       t_max),
      {}))...};
}

template<typename Series, int degree, int dimension, std::size_t... indices>
std::array<PoissonSeriesSubspace, sizeof...(indices)> AperiodicSeriesGenerator<
    Series, degree, dimension,
    std::index_sequence<indices...>>::Subspaces() {
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
      Instant const& t_min,
      Instant const& t_max);
  static std::array<PoissonSeriesSubspace, sizeof...(indices)> Subspaces(
      AngularFrequency const& ω);
};

template<typename Series, int degree, int dimension, std::size_t... indices>
std::array<Series, sizeof...(indices)> PeriodicSeriesGenerator<
    Series,
    degree, dimension,
    std::index_sequence<indices...>>::BasisElements(AngularFrequency const& ω,
                                                    Instant const& t_min,
                                                    Instant const& t_max) {
  Instant const t_mid = Barycentre({t_min, t_max});
  typename Series::AperiodicPolynomial const aperiodic_zero{{}, t_mid};
  return {Series(
      aperiodic_zero,
      {{ω,
        PolynomialGenerator<typename Series::PeriodicPolynomial, dimension>::
            template UnitPolynomials<typename Series::Polynomials, indices>(
                t_min, t_mid, t_max)}})...};
}

template<typename Series, int degree, int dimension, std::size_t... indices>
std::array<PoissonSeriesSubspace, sizeof...(indices)> PeriodicSeriesGenerator<
    Series, degree, dimension,
    std::index_sequence<indices...>>::Subspaces(AngularFrequency const& ω) {
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


template<typename Series, int degree>
auto PoissonSeriesBasisGenerator<Series, degree>::Basis(Instant const& t_min,
                                                        Instant const& t_max)
    -> std::array<Series, dimension * (degree + 1)> {
  return AperiodicSeriesGenerator<Series, degree, dimension>::BasisElements(
      t_min, t_max);
}

template<typename Series, int degree>
auto PoissonSeriesBasisGenerator<Series, degree>::Subspaces()
    -> std::array<PoissonSeriesSubspace, dimension*(degree + 1)> {
  return AperiodicSeriesGenerator<Series, degree, dimension>::Subspaces();
}

template<typename Series, int degree>
auto PoissonSeriesBasisGenerator<Series, degree>::Basis(
    AngularFrequency const& ω,
    Instant const& t_min,
    Instant const& t_max)
    -> std::array<Series, 2 * dimension*(degree + 1)> {
  return PeriodicSeriesGenerator<Series, degree, dimension>::BasisElements(
      ω, t_min, t_max);
}

template<typename Series, int degree>
auto PoissonSeriesBasisGenerator<Series, degree>::Subspaces(
    AngularFrequency const& ω)
    -> std::array<PoissonSeriesSubspace, 2 * dimension*(degree + 1)> {
  return PeriodicSeriesGenerator<Series, degree, dimension>::Subspaces(ω);
}

inline std::ostream& operator<<(std::ostream& out,
                                PoissonSeriesSubspace const& subspace) {
  return out << "{coordinate: " << static_cast<int>(subspace.coordinate_)
             << ", parity: " << static_cast<int>(subspace.parity_) << "}";
}

}  // namespace internal
}  // namespace _poisson_series_basis
}  // namespace numerics
}  // namespace principia
