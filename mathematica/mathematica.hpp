
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

class ExpressIn2 {
public:
  template<typename... Q>
  explicit ExpressIn2(Q const&... q) {
    FillUnits(q...);
  }

  template<typename Q>
  double operator()(Q const& q) const {
    return Divide<Q::Dimensions::Length>(
        Divide<Q::Dimensions::Mass>(
            Divide<Q::Dimensions::Time>(
                Divide<Q::Dimensions::Current>(
                    Divide<Q::Dimensions::Temperature>(
                        Divide<Q::Dimensions::Amount>(
                            Divide<Q::Dimensions::LuminousIntensity>(
                                Divide<Q::Dimensions::Angle>(
                                    q, std::get<std::optional<Angle>>(units_)),
                                std::get<std::optional<LuminousIntensity>>(
                                    units_)),
                            std::get<std::optional<Amount>>(units_)),
                        std::get<std::optional<Temperature>>(units_)),
                    std::get<std::optional<Current>>(units_)),
                std::get<std::optional<Time>>(units_)),
            std::get<std::optional<Mass>>(units_)),
        std::get<std::optional<Length>>(units_));
  }

 private:
  template<std::int64_t exponent, typename QLeft, typename QRight>
  static Quotient<QLeft, Exponentiation<QRight, exponent>> Divide(
    QLeft const& q_left,
    std::optional<QRight> const& q_right) {
    if constexpr (exponent == 0) {
      return q_left;
    }
    else {
      CHECK(q_right.has_value()) << "Missing unit to decompose " << q_left;
      return q_left / Pow<exponent>(q_right.value());
    }
  }


  template<typename QHead, typename... QTail>
  void FillUnits(QHead const& qhead, QTail const&... qtail) {
    // The following will fail to compile if QHead is not one of the units in
    // units_.
    auto& unit = std::get<std::optional<QHead>>(units_);
    CHECK(!unit.has_value()) << qhead << " is a redundant specification";
    unit = qhead;
    if constexpr (sizeof...(QTail) > 0) {
      FillUnits(qtail...);
    }
  }

  std::tuple<std::optional<Length>,
    std::optional<Mass>,
    std::optional<Time>,
    std::optional<Current>,
    std::optional<Temperature>,
    std::optional<Amount>,
    std::optional<LuminousIntensity>,
    std::optional<Angle>>
    units_;
};

std::string Apply(std::string const& function,
                  std::vector<std::string> const& arguments);

template<typename T>
std::string Option(std::string const& name, T const& right);

template<typename T>
std::string Assign(std::string const& name, T const& right);

template<typename T, typename U>
std::string PlottableDataset(std::vector<T> const& x, std::vector<U> const& y);

template<typename T>
std::string ToMathematica(std::vector<T> const& list);

template<typename It>
std::string ToMathematica(It begin, It end);

std::string ToMathematica(double const& real);

template<typename T, int size>
std::string ToMathematica(FixedVector<T, size> const& fixed_vector);

template<typename T>
std::string ToMathematica(R3Element<T> const& r3_element);

std::string ToMathematica(Quaternion const& quaternion);

template<typename D>
std::string ToMathematica(Quantity<D> const& quantity,
                          std::optional<ExpressIn2> express_in = std::nullopt);

template<typename S, typename F>
std::string ToMathematica(Vector<S, F> const& vector);

template<typename S, typename F>
std::string ToMathematica(Bivector<S, F> const& bivector);

template<typename V>
std::string ToMathematica(Point<V> const& point);

template<typename F>
std::string ToMathematica(DegreesOfFreedom<F> const& degrees_of_freedom);

template<typename... Types>
std::string ToMathematica(std::tuple<Types...> const& tuple);

template<typename R,
         typename = std::void_t<decltype(std::declval<R>().time)>,
         typename = std::void_t<decltype(std::declval<R>().degrees_of_freedom)>>
std::string ToMathematica(R ref);

std::string ToMathematica(
    astronomy::OrbitalElements::EquinoctialElements const& elements);

// Returns its argument.
std::string ToMathematica(std::string const& str);

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
