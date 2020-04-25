
#pragma once

#include <optional>
#include <string>
#include <tuple>
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
using quantities::Length;
using quantities::LuminousIntensity;
using quantities::Mass;
using quantities::Pow;
using quantities::Quantity;
using quantities::Quotient;
using quantities::Temperature;
using quantities::Time;

template<typename... Qs>
class ExpressIn2 {
 public:
  // static_assert

  static constexpr bool is_default = sizeof...(Qs) == 0;

  ExpressIn2(Qs const&... qs) : units_(std::make_tuple(qs...)) {}

  template<typename Q>
  double operator()(Q const& q) const {
    return Divide<Q::Dimensions::Length, Length>(
        Divide<Q::Dimensions::Mass, Mass>(
            Divide<Q::Dimensions::Time, Time>(
                Divide<Q::Dimensions::Current, Current>(
                    Divide<Q::Dimensions::Temperature, Temperature>(
                        Divide<Q::Dimensions::Amount, Amount>(
                            Divide<Q::Dimensions::LuminousIntensity,
                                   LuminousIntensity>(
                                Divide<Q::Dimensions::Angle, Angle>(q))))))));
  }

 //private:
  template<std::int64_t exponent, typename Q1, typename Q2>
  Quotient<Q2, Exponentiation<Q1, exponent>> Divide(Q2 const& q2) const {
    if constexpr (exponent == 0) {
      return q2;
    } else {
      return q2 / Pow<exponent>(std::get<Q1>(units_));
    }
  }

  std::tuple<Qs...> units_;
};

std::string Apply(std::string const& function,
                  std::vector<std::string> const& arguments);

template<typename T, typename... Qs>
std::string Option(std::string const& name,
                   T const& right,
                   ExpressIn2<Qs...> express_in = {});

template<typename T, typename... Qs>
std::string Assign(std::string const& name,
                   T const& right,
                   ExpressIn2<Qs...> express_in = {});

template<typename T, typename U, typename... Qs>
std::string PlottableDataset(std::vector<T> const& x,
                             std::vector<U> const& y,
                             ExpressIn2<Qs...> express_in = {});

template<typename T, typename... Qs>
std::string ToMathematica(std::vector<T> const& list,
                          ExpressIn2<Qs...> express_in = {});

template<typename It, typename... Qs>
std::string ToMathematica(It begin, It end,
                          ExpressIn2<Qs...> express_in = {});

template<typename... Qs>
std::string ToMathematica(double real,
                          ExpressIn2<Qs...> express_in = {});

template<typename T, int size, typename... Qs>
std::string ToMathematica(FixedVector<T, size> const& fixed_vector,
                          ExpressIn2<Qs...> express_in = {});

template<typename T, typename... Qs>
std::string ToMathematica(R3Element<T> const& r3_element,
                          ExpressIn2<Qs...> express_in = {});

template<typename... Qs>
std::string ToMathematica(Quaternion const& quaternion,
                          ExpressIn2<Qs...> express_in = {});

template<typename D, typename... Qs>
std::string ToMathematica(Quantity<D> const& quantity,
                          ExpressIn2<Qs...> express_in = {});

template<typename S, typename F, typename... Qs>
std::string ToMathematica(Vector<S, F> const& vector,
                          ExpressIn2<Qs...> express_in = {});

template<typename S, typename F, typename... Qs>
std::string ToMathematica(Bivector<S, F> const& bivector,
                          ExpressIn2<Qs...> express_in = {});

template<typename V, typename... Qs>
std::string ToMathematica(Point<V> const& point,
                          ExpressIn2<Qs...> express_in = {});

template<typename F, typename... Qs>
std::string ToMathematica(DegreesOfFreedom<F> const& degrees_of_freedom,
                          ExpressIn2<Qs...> express_in = {});

template<typename... Types>
std::string ToMathematica(std::tuple<Types...> const& tuple);

template<typename R,
         typename = std::void_t<decltype(std::declval<R>().time)>,
         typename = std::void_t<decltype(std::declval<R>().degrees_of_freedom)>>
std::string ToMathematica(R ref);

template<typename... Qs>
std::string ToMathematica(
    astronomy::OrbitalElements::EquinoctialElements const& elements,
    ExpressIn2<Qs...> express_in = {});

// Returns its argument.
template<typename... Qs>
std::string ToMathematica(std::string const& str,
                          ExpressIn2<Qs...> express_in = {});

// Wraps the string in quotes.
// TODO(egg): escape things properly.
std::string Escape(std::string const& str);

// TODO(phl): This doesn't work well for complex structures like orbital
// elements or trajectory iterators.  Surely we can do better.
template<typename T>
struct RemoveUnit;

template<typename T>
typename RemoveUnit<T>::Unitless ExpressIn(
    typename RemoveUnit<T>::Unit const& unit,
    T const& value);

template<typename V>
typename RemoveUnit<Point<V>>::Unitless ExpressIn(
    typename RemoveUnit<Point<V>>::Unit const& unit,
    Point<V> const& value);

template<typename T>
typename RemoveUnit<std::vector<T>>::Unitless ExpressIn(
    typename RemoveUnit<std::vector<T>>::Unit const& unit,
    std::vector<T> const& values);

}  // namespace internal_mathematica

using internal_mathematica::Apply;
using internal_mathematica::Assign;
using internal_mathematica::Escape;
using internal_mathematica::ExpressIn;
using internal_mathematica::ExpressIn2;
using internal_mathematica::Option;
using internal_mathematica::PlottableDataset;
using internal_mathematica::ToMathematica;

}  // namespace mathematica
}  // namespace principia

#include "mathematica/mathematica_body.hpp"
