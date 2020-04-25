
#pragma once

#include "mathematica/mathematica.hpp"

#include <cmath>
#include <string>
#include <tuple>
#include <vector>

#include "base/not_constructible.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace mathematica {
namespace internal_mathematica {

using astronomy::J2000;
using base::not_constructible;
using base::not_null;
using quantities::DebugString;
using quantities::IsFinite;
using quantities::si::Metre;
using quantities::si::Radian;
using quantities::si::Second;
namespace si = quantities::si;

// A helper struct to scan the elements of a tuple and stringify them.
template<int index, typename... Types>
struct TupleHelper : not_constructible {
  static void ToMathematicaStrings(std::tuple<Types...> const& tuple,
                                   std::vector<std::string>& expressions) {
    TupleHelper<index - 1, Types...>::ToMathematicaStrings(tuple, expressions);
    expressions.push_back(ToMathematica(std::get<index - 1>(tuple)));
  }
};

template<typename... Types>
struct TupleHelper<0, Types...> : not_constructible {
  static void ToMathematicaStrings(std::tuple<Types...> const& tuple,
                                   std::vector<std::string>& expressions) {}
};

inline std::string Apply(
    std::string const& function,
    std::vector<std::string> const& arguments) {
  std::string result;
  result += function;
  result += "[";
  for (int i = 0; i < arguments.size(); ++i) {
    result += arguments[i];
    result += (i + 1 < arguments.size() ? "," : "");
  }
  result += "]";
  return result;
}

template<typename T>
std::string Option(std::string const& name, T const& right) {
  return Apply("Rule", {name, ToMathematica(right)});
}

template<typename T>
std::string Assign(std::string const& name, T const& right) {
  return Apply("Set", {name, ToMathematica(right)}) + ";\n";
}

template<typename T, typename U>
std::string PlottableDataset(std::vector<T> const& x, std::vector<U> const& y) {
  std::vector<std::string> const xy = {ToMathematica(x), ToMathematica(y)};
  return Apply("Transpose", {ToMathematica(xy)});
}

template<typename T>
std::string ToMathematica(std::vector<T> const& list) {
  std::vector<std::string> expressions;
  expressions.reserve(list.size());
  for (auto const& expression : list) {
    expressions.emplace_back(ToMathematica(expression));
  }
  return Apply("List", expressions);
}

template<typename It>
std::string ToMathematica(It const begin, It const end) {
  std::vector<std::string> expressions;
  for (auto it = begin; it != end; ++it) {
    expressions.emplace_back(ToMathematica(*it));
  }
  return Apply("List", expressions);
}

inline std::string ToMathematica(double const& real) {
  if (std::isinf(real)) {
    if (real > 0.0) {
      return "Infinity";
    } else {
      return Apply("Minus", {"Infinity"});
    }
  } else if (std::isnan(real)) {
    return "Indeterminate";
  } else {
    std::string s = DebugString(real);
    s.replace(s.find("e"), 1, "*^");
    return Apply("SetPrecision", {s, "$MachinePrecision"});
  }
}

inline std::string ToMathematica(Quaternion const& quaternion) {
  return ToMathematica<double, R3Element<double>>(
      {quaternion.real_part(), quaternion.imaginary_part()});
}

template<typename T, int size>
std::string ToMathematica(FixedVector<T, size> const& fixed_vector) {
  std::vector<std::string> expressions;
  for (int i = 0; i < size; ++i) {
    expressions.emplace_back(ToMathematica(fixed_vector[i]));
  }
  return Apply("List", expressions);
}

template<typename T>
std::string ToMathematica(R3Element<T> const& r3_element) {
  std::vector<std::string> expressions;
  expressions.emplace_back(ToMathematica(r3_element.x));
  expressions.emplace_back(ToMathematica(r3_element.y));
  expressions.emplace_back(ToMathematica(r3_element.z));
  return Apply("List", expressions);
}

template<typename D, typename... Qs>
std::string ToMathematica(Quantity<D> const& quantity,
                          ExpressIn2<Qs...> express_in) {
  if (express_in.is_default()) {
    std::string s = DebugString(quantity);
    std::string const number = ToMathematica(quantity / si::Unit<Quantity<D>>);
    std::size_t const split = s.find(" ");
    std::string const units = Escape(s.substr(split, s.size()));
    return Apply("SetPrecision",
                 {Apply("Quantity", {number, units}), "$MachinePrecision"});
  } else {
    return ToMathematica(express_in(quantity));
  }
}

template<typename S, typename F>
std::string ToMathematica(Vector<S, F> const& vector) {
  return ToMathematica(vector.coordinates());
}

template<typename S, typename F>
std::string ToMathematica(Bivector<S, F> const& bivector) {
  return ToMathematica(bivector.coordinates());
}

template<typename V>
std::string ToMathematica(Point<V> const & point) {
  return ToMathematica(point - Point<V>());
}

template<typename F>
std::string ToMathematica(DegreesOfFreedom<F> const& degrees_of_freedom) {
  return Apply(
      "List",
      std::vector<std::string>{
          ToMathematica(ExpressIn(Metre, degrees_of_freedom.position())),
          ToMathematica(
              ExpressIn(Metre / Second, degrees_of_freedom.velocity()))});
}

template<typename... Types>
std::string ToMathematica(std::tuple<Types...> const& tuple) {
  std::vector<std::string> expressions;
  expressions.reserve(sizeof...(Types));
  TupleHelper<sizeof...(Types), Types...>::ToMathematicaStrings(
      tuple, expressions);
  return Apply("List", expressions);
}

template<typename R, typename, typename>
std::string ToMathematica(R const ref) {
  return Apply(
      "List",
      std::vector<std::string>{ToMathematica(ExpressIn(Second, ref.time)),
                               ToMathematica(ref.degrees_of_freedom)});
}

inline std::string ToMathematica(
    astronomy::OrbitalElements::EquinoctialElements const& elements) {
  return ToMathematica(std::make_tuple((elements.t - J2000) / Second,
                                       elements.a / Metre,
                                       elements.h,
                                       elements.k,
                                       elements.λ / Radian,
                                       elements.p,
                                       elements.q,
                                       elements.pʹ,
                                       elements.qʹ));
}

inline std::string ToMathematica(std::string const& str) {
  return str;
}

// Wraps the string in quotes.
// TODO(egg): escape things properly.
inline std::string Escape(std::string const& str) {
  std::string result = {"\""};
  result += str;
  result += "\"";
  return result;
}

template<typename T>
struct RemoveUnit<Quantity<T>> {
  using Unit = Quantity<T>;
  using Unitless = double;
};

template<typename T, typename F>
struct RemoveUnit<Vector<T, F>> {
  using Unit = typename RemoveUnit<T>::Unit;
  using Unitless = Vector<typename RemoveUnit<T>::Unitless, F>;
};

template<typename V>
struct RemoveUnit<Point<V>> {
  using Unit = typename RemoveUnit<V>::Unit;
  using Unitless = Point<typename RemoveUnit<V>::Unitless>;
};

template<typename T>
struct RemoveUnit<std::vector<T>> {
  using Unit = typename RemoveUnit<T>::Unit;
  using Unitless = std::vector<typename RemoveUnit<T>::Unitless>;
};

template<typename T>
typename RemoveUnit<T>::Unitless ExpressIn(
    typename RemoveUnit<T>::Unit const& unit,
    T const& value) {
  return value / unit;
}

template<typename V>
typename RemoveUnit<Point<V>>::Unitless ExpressIn(
    typename RemoveUnit<Point<V>>::Unit const& unit,
    Point<V> const& value) {
  return (value - Point<V>{}) / unit +
         typename RemoveUnit<Point<V>>::Unitless{};
}

template<typename T>
typename RemoveUnit<std::vector<T>>::Unitless ExpressIn(
    typename RemoveUnit<std::vector<T>>::Unit const& unit,
    std::vector<T> const& values) {
  typename RemoveUnit<std::vector<T>>::Unitless result;
  for (auto const& value : values) {
    result.push_back(ExpressIn(unit, value));
  }
  return result;
}

}  // namespace internal_mathematica
}  // namespace mathematica
}  // namespace principia
