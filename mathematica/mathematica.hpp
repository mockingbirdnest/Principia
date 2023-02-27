#pragma once

#include <map>
#include <string>
#include <tuple>
#include <type_traits>
#include <vector>

#include "astronomy/orbital_elements.hpp"
#include "base/traits.hpp"
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
#include "numerics/polynomial.hpp"
#include "numerics/unbounded_arrays.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "quantities/elementary_functions.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "quantities/tuples.hpp"

namespace principia {
namespace mathematica {
namespace internal_mathematica {

using numerics::DoublePrecision;
using numerics::FixedVector;
using numerics::PiecewisePoissonSeries;
using numerics::PoissonSeries;
using numerics::PolynomialInMonomialBasis;
using numerics::UnboundedLowerTriangularMatrix;
using numerics::UnboundedUpperTriangularMatrix;
using numerics::UnboundedVector;
using physics::DegreesOfFreedom;
using quantities::Amount;
using quantities::Angle;
using quantities::Current;
using quantities::Exponentiation;
using quantities::is_tuple_v;
using quantities::Length;
using quantities::LuminousIntensity;
using quantities::Mass;
using quantities::Pow;
using quantities::Quantity;
using quantities::Quotient;
using quantities::Temperature;
using quantities::Time;
using namespace principia::base::_traits;
using namespace principia::geometry::_grassmann;
using namespace principia::geometry::_point;
using namespace principia::geometry::_quaternion;
using namespace principia::geometry::_r3_element;
using namespace principia::geometry::_r3x3_matrix;
using namespace principia::geometry::_symmetric_bilinear_form;
namespace si = quantities::si;

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
                          OptionalExpressIn express_in = std::nullopt);

template<typename T, int size, typename OptionalExpressIn = std::nullopt_t>
std::string ToMathematica(FixedVector<T, size> const& fixed_vector,
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

template<typename Tuple,
         typename = std::enable_if_t<is_tuple_v<Tuple>>,
         typename OptionalExpressIn = std::nullopt_t>
std::string ToMathematica(Tuple const& tuple,
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

template<typename R,
         typename = std::void_t<decltype(std::declval<R>().time)>,
         typename = std::void_t<decltype(std::declval<R>().degrees_of_freedom)>,
         typename OptionalExpressIn = std::nullopt_t>
std::string ToMathematica(R ref,
                          OptionalExpressIn express_in = std::nullopt);

template<typename V, typename A, int d,
         template<typename, typename, int> class E,
         typename OptionalExpressIn = std::nullopt_t>
std::string ToMathematicaBody(
    PolynomialInMonomialBasis<V, A, d, E> const& polynomial,
    OptionalExpressIn express_in);

template<typename V, typename A, int d,
         template<typename, typename, int> class E,
         typename OptionalExpressIn = std::nullopt_t>
std::string ToMathematica(
    PolynomialInMonomialBasis<V, A, d, E> const& polynomial,
    OptionalExpressIn express_in = std::nullopt);


template<typename V, int ad, int pd,
         template<typename, typename, int> class E,
         typename OptionalExpressIn = std::nullopt_t>
std::string ToMathematicaBody(
    PoissonSeries<V, ad, pd, E> const& series,
    OptionalExpressIn express_in);

template<typename V, int ad, int pd,
         template<typename, typename, int> class E,
         typename OptionalExpressIn = std::nullopt_t>
std::string ToMathematica(PoissonSeries<V, ad, pd, E> const& series,
                          OptionalExpressIn express_in = std::nullopt);

template<typename V, int ad, int pd,
         template<typename, typename, int> class E,
         typename OptionalExpressIn = std::nullopt_t>
std::string ToMathematicaBody(
    PiecewisePoissonSeries<V, ad, pd, E> const& series,
    OptionalExpressIn express_in);

template<typename V, int ad, int pd,
         template<typename, typename, int> class E,
         typename OptionalExpressIn = std::nullopt_t>
std::string ToMathematica(PiecewisePoissonSeries<V, ad, pd, E> const& series,
                          OptionalExpressIn express_in = std::nullopt);

template<typename OptionalExpressIn = std::nullopt_t>
std::string ToMathematica(
    astronomy::OrbitalElements::EquinoctialElements const& elements,
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

}  // namespace internal_mathematica

using internal_mathematica::Apply;
using internal_mathematica::Evaluate;
using internal_mathematica::ExpressIn;
using internal_mathematica::ExpressInSIUnits;
using internal_mathematica::PlottableDataset;
using internal_mathematica::PreserveUnits;
using internal_mathematica::Rule;
using internal_mathematica::Set;
using internal_mathematica::ToMathematica;
using internal_mathematica::ToMathematicaBody;

}  // namespace mathematica
}  // namespace principia

#include "mathematica/mathematica_body.hpp"
