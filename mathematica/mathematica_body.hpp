
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

template<typename T, typename... Qs>
std::string Option(std::string const& name,
                   T const& right,
                   ExpressIn2<Qs...> express_in) {
  return Apply("Rule", {name, ToMathematica(right, express_in)});
}

template<typename T, typename... Qs>
std::string Assign(std::string const& name,
                   T const& right,
                   ExpressIn2<Qs...> express_in) {
  return Apply("Set", {name, ToMathematica(right, express_in)}) + ";\n";
}

template<typename T, typename U, typename... Qs>
std::string PlottableDataset(std::vector<T> const& x,
                             std::vector<U> const& y,
                             ExpressIn2<Qs...> express_in) {
  std::vector<std::string> const xy = {ToMathematica(x, express_in),
                                       ToMathematica(y, express_in)};
  return Apply("Transpose", {ToMathematica(xy, express_in)});
}

template<typename T, typename... Qs>
std::string ToMathematica(std::vector<T> const& list,
                          ExpressIn2<Qs...> express_in) {
  std::vector<std::string> expressions;
  expressions.reserve(list.size());
  for (auto const& expression : list) {
    expressions.emplace_back(ToMathematica(expression, express_in));
  }
  return Apply("List", expressions);
}

template<typename It, typename... Qs>
std::string ToMathematica(It const begin, It const end,
                          ExpressIn2<Qs...> express_in) {
  std::vector<std::string> expressions;
  for (auto it = begin; it != end; ++it) {
    expressions.emplace_back(ToMathematica(*it, express_in));
  }
  return Apply("List", expressions);
}

template<typename... Qs>
std::string ToMathematica(double const real,
                          ExpressIn2<Qs...> /*express_in*/) {
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

template<typename... Qs>
std::string ToMathematica(Quaternion const& quaternion,
                          ExpressIn2<Qs...> express_in) {
  return ToMathematica(
      std::tuple{quaternion.real_part(), quaternion.imaginary_part()});
}

template<typename T, int size, typename... Qs>
std::string ToMathematica(FixedVector<T, size> const& fixed_vector,
                          ExpressIn2<Qs...> express_in) {
  std::vector<std::string> expressions;
  for (int i = 0; i < size; ++i) {
    expressions.emplace_back(ToMathematica(fixed_vector[i], express_in));
  }
  return Apply("List", expressions);
}

template<typename T, typename... Qs>
std::string ToMathematica(R3Element<T> const& r3_element,
                          ExpressIn2<Qs...> express_in) {
  std::vector<std::string> expressions;
  expressions.emplace_back(ToMathematica(r3_element.x, express_in));
  expressions.emplace_back(ToMathematica(r3_element.y, express_in));
  expressions.emplace_back(ToMathematica(r3_element.z, express_in));
  return Apply("List", expressions);
}

template<typename D, typename... Qs>
std::string ToMathematica(Quantity<D> const& quantity,
                          ExpressIn2<Qs...> express_in) {
  if constexpr (ExpressIn2<Qs...>::is_default) {
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

template<typename S, typename F, typename... Qs>
std::string ToMathematica(Vector<S, F> const& vector,
                          ExpressIn2<Qs...> express_in) {
  return ToMathematica(vector.coordinates(), express_in);
}

template<typename S, typename F, typename... Qs>
std::string ToMathematica(Bivector<S, F> const& bivector,
                          ExpressIn2<Qs...> express_in) {
  return ToMathematica(bivector.coordinates(), express_in);
}

template<typename V, typename... Qs>
std::string ToMathematica(Point<V> const & point,
                          ExpressIn2<Qs...> express_in) {
  return ToMathematica(point - Point<V>(), express_in);
}

template<typename F, typename... Qs>
std::string ToMathematica(DegreesOfFreedom<F> const& degrees_of_freedom,
                          ExpressIn2<Qs...> express_in) {
  return Apply(
      "List",
      std::vector<std::string>{
          ToMathematica(degrees_of_freedom.position(), express_in),
          ToMathematica(degrees_of_freedom.velocity(), express_in)});
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

template<typename... Qs>
std::string ToMathematica(
    astronomy::OrbitalElements::EquinoctialElements const& elements,
    ExpressIn2<Qs...> express_in) {
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

template<typename... Qs>
std::string ToMathematica(std::string const& str,
                          ExpressIn2<Qs...>) {
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
