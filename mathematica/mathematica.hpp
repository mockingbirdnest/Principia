
#pragma once

#include <string>
#include <tuple>
#include <type_traits>
#include <vector>

#include "astronomy/orbital_elements.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/point.hpp"
#include "geometry/quaternion.hpp"
#include "geometry/r3_element.hpp"
#include "numerics/fixed_arrays.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "quantities/elementary_functions.hpp"
#include "quantities/quantities.hpp"
#include "quantities/tuples.hpp"

namespace principia {
namespace mathematica {
namespace internal_mathematica {

using geometry::Bivector;
using geometry::Point;
using geometry::Quaternion;
using geometry::R3Element;
using geometry::Vector;
using numerics::FixedVector;
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

// A helper class for type erasure of quantities.  It may be used with the
// functions in this file to remove the dimensions of quantities (we know that
// Mathematica is sluggish when processing quantities).  Usage:
//
//   ToMathematica(... , ExpressIn(Metre, Second, Degree));
//
// The construction parameters must be for the base SI units.  They may be in
// any order.  If not enough parameters are specified for complete type erasure
// of the argument, compilation fails.  An object with no template parameters
// may be used to indicate that type erasure should not happen.
template<typename... Qs>
class ExpressIn {
 public:
  static constexpr bool is_default = sizeof...(Qs) == 0;

  // Check that only base SI units are specified.  If a unit is specified
  // multiple times, the calls to std::get in operator() won't compile.
  static_assert(
      is_default ||
      ((std::is_same_v<Qs, Length> || std::is_same_v<Qs, Mass> ||
        std::is_same_v<Qs, Time> || std::is_same_v<Qs, Current> ||
        std::is_same_v<Qs, Temperature> || std::is_same_v<Qs, Amount> ||
        std::is_same_v<Qs, LuminousIntensity> || std::is_same_v<Qs, Angle>) ||
       ...));

  ExpressIn(Qs const&... qs);

  template<typename Q>
  double operator()(Q const& q) const;

 private:
  template<std::int64_t exponent, typename Q1, typename Q2>
  Quotient<Q2, Exponentiation<Q1, exponent>> Divide(Q2 const& q2) const;

  std::tuple<Qs...> units_;
};

std::string Apply(std::string const& function,
                  std::vector<std::string> const& arguments);

template<typename T, typename... Qs>
std::string Option(std::string const& name,
                   T const& right,
                   ExpressIn<Qs...> express_in = {});

template<typename T, typename... Qs>
std::string Assign(std::string const& name,
                   T const& right,
                   ExpressIn<Qs...> express_in = {});

template<typename T, typename U, typename... Qs>
std::string PlottableDataset(std::vector<T> const& x,
                             std::vector<U> const& y,
                             ExpressIn<Qs...> express_in = {});

template<typename T, typename... Qs>
std::string ToMathematica(std::vector<T> const& list,
                          ExpressIn<Qs...> express_in = {});

template<typename It, typename... Qs>
std::string ToMathematica(It begin, It end,
                          ExpressIn<Qs...> express_in = {});

template<typename... Qs>
std::string ToMathematica(double real,
                          ExpressIn<Qs...> express_in = {});

template<typename T, int size, typename... Qs>
std::string ToMathematica(FixedVector<T, size> const& fixed_vector,
                          ExpressIn<Qs...> express_in = {});

template<typename T, typename... Qs>
std::string ToMathematica(R3Element<T> const& r3_element,
                          ExpressIn<Qs...> express_in = {});

template<typename... Qs>
std::string ToMathematica(Quaternion const& quaternion,
                          ExpressIn<Qs...> express_in = {});

template<typename D, typename... Qs>
std::string ToMathematica(Quantity<D> const& quantity,
                          ExpressIn<Qs...> express_in = {});

template<typename S, typename F, typename... Qs>
std::string ToMathematica(Vector<S, F> const& vector,
                          ExpressIn<Qs...> express_in = {});

template<typename S, typename F, typename... Qs>
std::string ToMathematica(Bivector<S, F> const& bivector,
                          ExpressIn<Qs...> express_in = {});

template<typename V, typename... Qs>
std::string ToMathematica(Point<V> const& point,
                          ExpressIn<Qs...> express_in = {});

template<typename F, typename... Qs>
std::string ToMathematica(DegreesOfFreedom<F> const& degrees_of_freedom,
                          ExpressIn<Qs...> express_in = {});

template<typename Tuple,
         typename = std::enable_if_t<is_tuple_v<Tuple>>,
         typename... Qs>
std::string ToMathematica(Tuple const& tuple,
                          ExpressIn<Qs...> express_in = {});

template<typename R,
         typename = std::void_t<decltype(std::declval<R>().time)>,
         typename = std::void_t<decltype(std::declval<R>().degrees_of_freedom)>,
         typename... Qs>
std::string ToMathematica(R ref,
                          ExpressIn<Qs...> express_in = {});

template<typename... Qs>
std::string ToMathematica(
    astronomy::OrbitalElements::EquinoctialElements const& elements,
    ExpressIn<Qs...> express_in = {});

// Returns its argument.
template<typename... Qs>
std::string ToMathematica(std::string const& str,
                          ExpressIn<Qs...> express_in = {});

// Wraps the string in quotes.
// TODO(egg): escape things properly.
std::string Escape(std::string const& str);

}  // namespace internal_mathematica

using internal_mathematica::Apply;
using internal_mathematica::Assign;
using internal_mathematica::Escape;
using internal_mathematica::ExpressIn;
using internal_mathematica::Option;
using internal_mathematica::PlottableDataset;
using internal_mathematica::ToMathematica;

}  // namespace mathematica
}  // namespace principia

#include "mathematica/mathematica_body.hpp"
