#pragma once

#include <map>
#include <string>
#include <tuple>
#include <type_traits>
#include <vector>

#include "astronomy/orbital_elements.hpp"
#include "base/traits.hpp"
#include "boost/multiprecision/cpp_bin_float.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/point.hpp"
#include "geometry/quaternion.hpp"
#include "geometry/r3_element.hpp"
#include "geometry/r3x3_matrix.hpp"
#include "geometry/symmetric_bilinear_form.hpp"
#include "numerics/double_precision.hpp"
#include "numerics/fixed_arrays.hpp"
#include "numerics/piecewise_poisson_series.hpp"
#include "numerics/poisson_series.hpp"
#include "numerics/polynomial_in_monomial_basis.hpp"
#include "numerics/polynomial_in_чебышёв_basis.hpp"
#include "numerics/unbounded_arrays.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "physics/discrete_trajectory.hpp"
#include "quantities/elementary_functions.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "quantities/tuples.hpp"

namespace principia {
namespace mathematica {
namespace _mathematica {
namespace internal {

using namespace boost::multiprecision;
using namespace principia::astronomy::_orbital_elements;
using namespace principia::base::_traits;
using namespace principia::geometry::_grassmann;
using namespace principia::geometry::_point;
using namespace principia::geometry::_quaternion;
using namespace principia::geometry::_r3_element;
using namespace principia::geometry::_r3x3_matrix;
using namespace principia::geometry::_symmetric_bilinear_form;
using namespace principia::numerics::_double_precision;
using namespace principia::numerics::_fixed_arrays;
using namespace principia::numerics::_piecewise_poisson_series;
using namespace principia::numerics::_poisson_series;
using namespace principia::numerics::_polynomial_in_monomial_basis;
using namespace principia::numerics::_polynomial_in_чебышёв_basis;
using namespace principia::numerics::_unbounded_arrays;
using namespace principia::physics::_degrees_of_freedom;
using namespace principia::physics::_discrete_trajectory;
using namespace principia::quantities::_elementary_functions;
using namespace principia::quantities::_named_quantities;
using namespace principia::quantities::_quantities;
using namespace principia::quantities::_si;
using namespace principia::quantities::_tuples;

template<typename F>
using DiscreteTrajectoryValueType =
    principia::physics::_discrete_trajectory_types::internal::value_type<F>;

// A helper class for type erasure of quantities.  It may be used with the
// functions in this file to remove the dimensions of quantities (we know that
// Mathematica is sluggish when processing quantities).  Usage:
//
//   ToMathematica(... , ExpressIn(Metre, Second, Degree));
//
// The construction parameters must be values of distinct SI base quantities
// (or Angle). They define a system of units.  They may be in any order.  If the
// other arguments of the functions contain quantities that are not spanned by
// that system of units, the call is ill-formed.
//
// A shortcut is provided for SI units:
//
//   ToMathematica(..., ExpressInSIUnits);
//
// and for preserving the units in the Mathematica output:
//
//   ToMathematica(..., PreserveUnits);
//
// One of ExpressIn... or PreserveUnits must be present as soon as the data
// being converted contains (dimensionful) quantities.

template<typename... Qs>
class ExpressIn {
 public:
  // Check that only SI base quantities or Angle are specified.
  static_assert(
      ((std::is_same_v<Qs, Length> || std::is_same_v<Qs, Mass> ||
        std::is_same_v<Qs, Time> || std::is_same_v<Qs, Current> ||
        std::is_same_v<Qs, Temperature> || std::is_same_v<Qs, Amount> ||
        std::is_same_v<Qs, LuminousIntensity> || std::is_same_v<Qs, Angle>) &&
       ...), "Must instantiate with SI base quantities or Angle");

  // Check that all quantities are different.
  static_assert(all_different_v<Qs...>, "Must use different quantities");

  constexpr ExpressIn(Qs const&... qs)  // NOLINT(runtime/explicit)
      : units_(std::make_tuple(qs...)) {}

  template<typename Q>
  double operator()(Q const& q) const;

 private:
  template<std::int64_t exponent, typename Q1, typename Q2>
  Quotient<Q2, Exponentiation<Q1, exponent>> Divide(Q2 const& q2) const;

  std::tuple<Qs...> units_;
};

constexpr auto ExpressInSIUnits = ExpressIn(si::Unit<Length>,
                                            si::Unit<Mass>,
                                            si::Unit<Time>,
                                            si::Unit<Current>,
                                            si::Unit<Temperature>,
                                            si::Unit<Amount>,
                                            si::Unit<LuminousIntensity>,
                                            si::Unit<Angle>);
constexpr struct {} PreserveUnits;

template<typename T, typename OptionalExpressIn = std::nullopt_t>
std::string Apply(std::string const& name,
                  T const& right,
                  OptionalExpressIn express_in = std::nullopt);

std::string Evaluate(std::string const& expression);

template<typename T, typename OptionalExpressIn = std::nullopt_t>
std::string Rule(std::string const& name,
                 T const& right,
                 OptionalExpressIn express_in = std::nullopt);

template<typename T, typename OptionalExpressIn = std::nullopt_t>
std::string Set(std::string const& name,
                T const& right,
                OptionalExpressIn express_in = std::nullopt);

template<typename T, typename U, typename OptionalExpressIn = std::nullopt_t>
std::string PlottableDataset(std::vector<T> const& x,
                             std::vector<U> const& y,
                             OptionalExpressIn express_in = std::nullopt);

template<typename T, typename OptionalExpressIn = std::nullopt_t>
std::string ToMathematica(std::vector<T> const& list,
                          OptionalExpressIn express_in = std::nullopt);

template<typename It, typename OptionalExpressIn = std::nullopt_t>
std::string ToMathematica(It begin, It end,
                          OptionalExpressIn express_in = std::nullopt);

template<typename OptionalExpressIn = std::nullopt_t>
std::string ToMathematica(bool b,
                          OptionalExpressIn express_in = std::nullopt);

template<typename T,
         typename = std::enable_if_t<std::is_integral_v<T>>,
         typename OptionalExpressIn = std::nullopt_t>
std::string ToMathematica(T integer,
                          OptionalExpressIn express_in = std::nullopt);

template<typename T,
         typename = std::enable_if_t<std::is_floating_point_v<T>>,
         typename OptionalExpressIn = std::nullopt_t,
         typename = void>
std::string ToMathematica(T real,
                          OptionalExpressIn express_in = std::nullopt,
                          std::int64_t base = 16);

template<unsigned digits,
         typename OptionalExpressIn = std::nullopt_t>
std::string ToMathematica(
    number<backends::cpp_bin_float<digits>> const& cpp_bin_float,
    OptionalExpressIn express_in = std::nullopt,
    std::int64_t base = 16);

template<typename T, std::int64_t size,
         typename OptionalExpressIn = std::nullopt_t>
std::string ToMathematica(FixedVector<T, size> const& fixed_vector,
                          OptionalExpressIn express_in = std::nullopt);

template<typename T, std::int64_t rows, std::int64_t columns,
         typename OptionalExpressIn = std::nullopt_t>
std::string ToMathematica(FixedMatrix<T, rows, columns> const& fixed_matrix,
                          OptionalExpressIn express_in = std::nullopt);

template<typename T, typename OptionalExpressIn = std::nullopt_t>
std::string ToMathematica(R3Element<T> const& r3_element,
                          OptionalExpressIn express_in = std::nullopt);

template<typename T, typename OptionalExpressIn = std::nullopt_t>
std::string ToMathematica(R3x3Matrix<T> const& r3x3_matrix,
                          OptionalExpressIn express_in = std::nullopt);

template<typename OptionalExpressIn = std::nullopt_t>
std::string ToMathematica(Quaternion const& quaternion,
                          OptionalExpressIn express_in = std::nullopt);

template<typename D, typename... Qs>
std::string ToMathematica(Quantity<D> const& quantity,
                          ExpressIn<Qs...> express_in);

template<typename D>
std::string ToMathematica(Quantity<D> const& quantity,
                          decltype(PreserveUnits) express_in);

template<typename D>
std::string ToMathematica(Quantity<D> const& quantity,
                          std::nullopt_t express_in = std::nullopt);

template<typename T, typename OptionalExpressIn = std::nullopt_t>
std::string ToMathematica(DoublePrecision<T> const& double_precision,
                          OptionalExpressIn express_in = std::nullopt);

template<typename S, typename F, typename OptionalExpressIn = std::nullopt_t>
std::string ToMathematica(Vector<S, F> const& vector,
                          OptionalExpressIn express_in = std::nullopt);

template<typename S, typename F, typename OptionalExpressIn = std::nullopt_t>
std::string ToMathematica(Bivector<S, F> const& bivector,
                          OptionalExpressIn express_in = std::nullopt);

template<typename V, typename OptionalExpressIn = std::nullopt_t>
std::string ToMathematica(Point<V> const& point,
                          OptionalExpressIn express_in = std::nullopt);

template<typename S, typename F,
         template<typename, typename> typename M,
         typename OptionalExpressIn = std::nullopt_t>
std::string ToMathematica(SymmetricBilinearForm<S, F, M> const& form,
                          OptionalExpressIn express_in = std::nullopt);

template<typename F, typename OptionalExpressIn = std::nullopt_t>
std::string ToMathematica(DegreesOfFreedom<F> const& degrees_of_freedom,
                          OptionalExpressIn express_in = std::nullopt);

template<typename F, typename OptionalExpressIn = std::nullopt_t>
std::string ToMathematica(
    RelativeDegreesOfFreedom<F> const& relative_degrees_of_freedom,
    OptionalExpressIn express_in = std::nullopt);

template<typename Tuple,
         typename = std::enable_if_t<is_tuple_v<Tuple>>,
         typename OptionalExpressIn = std::nullopt_t>
std::string ToMathematica(Tuple const& tuple,
                          OptionalExpressIn express_in = std::nullopt);

template<typename Scalar, typename OptionalExpressIn = std::nullopt_t>
std::string ToMathematica(UnboundedMatrix<Scalar> const& matrix,
                          OptionalExpressIn express_in = std::nullopt);

template<typename Scalar, typename OptionalExpressIn = std::nullopt_t>
std::string ToMathematica(UnboundedLowerTriangularMatrix<Scalar> const& matrix,
                          OptionalExpressIn express_in = std::nullopt);

template<typename Scalar, typename OptionalExpressIn = std::nullopt_t>
std::string ToMathematica(UnboundedUpperTriangularMatrix<Scalar> const& matrix,
                          OptionalExpressIn express_in = std::nullopt);

template<typename Scalar, typename OptionalExpressIn = std::nullopt_t>
std::string ToMathematica(UnboundedVector<Scalar> const& vector,
                          OptionalExpressIn express_in = std::nullopt);

template<typename F,
         typename OptionalExpressIn = std::nullopt_t>
std::string ToMathematica(DiscreteTrajectory<F> const& trajectory,
                          OptionalExpressIn express_in = std::nullopt);

template<typename F,
         typename OptionalExpressIn = std::nullopt_t>
std::string ToMathematica(DiscreteTrajectoryValueType<F> const& v,
                          OptionalExpressIn express_in = std::nullopt);

template<typename V, typename A, int d,
         typename OptionalExpressIn = std::nullopt_t>
std::string ToMathematicaBody(
    PolynomialInMonomialBasis<V, A, d> const& polynomial,
    OptionalExpressIn express_in);

template<typename V, typename A, int d,
         typename OptionalExpressIn = std::nullopt_t>
std::string ToMathematica(
    PolynomialInMonomialBasis<V, A, d> const& polynomial,
    OptionalExpressIn express_in = std::nullopt);

template<typename V, typename A, int d,
         typename OptionalExpressIn = std::nullopt_t>
std::string ToMathematicaBody(
    PolynomialInЧебышёвBasis<V, A, d> const& polynomial,
    OptionalExpressIn express_in);

template<typename V, typename A, int d,
         typename OptionalExpressIn = std::nullopt_t>
std::string ToMathematica(
    PolynomialInЧебышёвBasis<V, A, d> const& polynomial,
    OptionalExpressIn express_in = std::nullopt);

template<typename V, int ad, int pd,
         typename OptionalExpressIn = std::nullopt_t>
std::string ToMathematicaBody(
    PoissonSeries<V, ad, pd> const& series,
    OptionalExpressIn express_in);

template<typename V, int ad, int pd,
         typename OptionalExpressIn = std::nullopt_t>
std::string ToMathematica(PoissonSeries<V, ad, pd> const& series,
                          OptionalExpressIn express_in = std::nullopt);

template<typename V, int ad, int pd,
         typename OptionalExpressIn = std::nullopt_t>
std::string ToMathematicaBody(
    PiecewisePoissonSeries<V, ad, pd> const& series,
    OptionalExpressIn express_in);

template<typename V, int ad, int pd,
         typename OptionalExpressIn = std::nullopt_t>
std::string ToMathematica(PiecewisePoissonSeries<V, ad, pd> const& series,
                          OptionalExpressIn express_in = std::nullopt);

template<typename OptionalExpressIn = std::nullopt_t>
std::string ToMathematica(OrbitalElements::EquinoctialElements const& elements,
                          OptionalExpressIn express_in = std::nullopt);

template<typename T, typename OptionalExpressIn = std::nullopt_t>
std::string ToMathematica(std::optional<T> const& opt,
                          OptionalExpressIn express_in = std::nullopt);

template<typename OptionalExpressIn = std::nullopt_t>
std::string ToMathematica(char const* str,
                          OptionalExpressIn express_in = std::nullopt);
template<typename OptionalExpressIn = std::nullopt_t>
std::string ToMathematica(std::string const& str,
                          OptionalExpressIn express_in = std::nullopt);

// Does not wrap its arguments in ToMathematica.  Do not export.  Do not use
// externally.
std::string RawApply(std::string const& function,
                     std::vector<std::string> const& arguments);

}  // namespace internal

using internal::Apply;
using internal::Evaluate;
using internal::ExpressIn;
using internal::ExpressInSIUnits;
using internal::PlottableDataset;
using internal::PreserveUnits;
using internal::Rule;
using internal::Set;
using internal::ToMathematica;
using internal::ToMathematicaBody;

}  // namespace _mathematica
}  // namespace mathematica
}  // namespace principia

#include "mathematica/mathematica_body.hpp"
