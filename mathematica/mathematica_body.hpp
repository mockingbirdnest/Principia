
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
template<int index, typename Tuple, typename... Qs>
struct TupleHelper : not_constructible {
  static void ToMathematicaStrings(Tuple const& tuple,
                                   std::vector<std::string>& expressions,
                                   ExpressIn<Qs...> express_in) {
    TupleHelper<index - 1, Tuple, Qs...>::ToMathematicaStrings(
        tuple, expressions, express_in);
    expressions.push_back(ToMathematica(std::get<index - 1>(tuple),
                          express_in));
  }
};

template<typename Tuple, typename... Qs>
struct TupleHelper<0, Tuple, Qs...> : not_constructible {
  static void ToMathematicaStrings(Tuple const& tuple,
                                   std::vector<std::string>& expressions,
                                   ExpressIn<Qs...> express_in) {}
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

template<typename... Qs>
ExpressIn<Qs...>::ExpressIn(Qs const&... qs)
    : units_(std::make_tuple(qs...)) {}

template<typename... Qs>
template<typename Q>
double ExpressIn<Qs...>::operator()(Q const& q) const {
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

template<typename... Qs>
template<std::int64_t exponent, typename Q1, typename Q2>
Quotient<Q2, Exponentiation<Q1, exponent>> ExpressIn<Qs...>::Divide(
    Q2 const& q2) const {
  if constexpr (exponent == 0) {
    return q2;
  } else {
    return q2 / Pow<exponent>(std::get<Q1>(units_));
  }
}

template<typename T, typename... Qs>
std::string Option(std::string const& name,
                   T const& right,
                   ExpressIn<Qs...> express_in) {
  return Apply("Rule", {name, ToMathematica(right, express_in)});
}

template<typename T, typename... Qs>
std::string Assign(std::string const& name,
                   T const& right,
                   ExpressIn<Qs...> express_in) {
  return Apply("Set", {name, ToMathematica(right, express_in)}) + ";\n";
}

template<typename T, typename U, typename... Qs>
std::string PlottableDataset(std::vector<T> const& x,
                             std::vector<U> const& y,
                             ExpressIn<Qs...> express_in) {
  std::vector<std::string> const xy = {ToMathematica(x, express_in),
                                       ToMathematica(y, express_in)};
  return Apply("Transpose", {ToMathematica(xy, express_in)});
}

template<typename T, typename... Qs>
std::string ToMathematica(std::vector<T> const& list,
                          ExpressIn<Qs...> express_in) {
  std::vector<std::string> expressions;
  expressions.reserve(list.size());
  for (auto const& expression : list) {
    expressions.emplace_back(ToMathematica(expression, express_in));
  }
  return Apply("List", expressions);
}

template<typename It, typename... Qs>
std::string ToMathematica(It const begin, It const end,
                          ExpressIn<Qs...> express_in) {
  std::vector<std::string> expressions;
  for (auto it = begin; it != end; ++it) {
    expressions.emplace_back(ToMathematica(*it, express_in));
  }
  return Apply("List", expressions);
}

template<typename... Qs>
std::string ToMathematica(double const real,
                          ExpressIn<Qs...> /*express_in*/) {
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
                          ExpressIn<Qs...> express_in) {
  return ToMathematica(
      std::tuple{quaternion.real_part(), quaternion.imaginary_part()},
      express_in);
}

template<typename T, int size, typename... Qs>
std::string ToMathematica(FixedVector<T, size> const& fixed_vector,
                          ExpressIn<Qs...> express_in) {
  std::vector<std::string> expressions;
  for (int i = 0; i < size; ++i) {
    expressions.emplace_back(ToMathematica(fixed_vector[i], express_in));
  }
  return Apply("List", expressions);
}

template<typename T, typename... Qs>
std::string ToMathematica(R3Element<T> const& r3_element,
                          ExpressIn<Qs...> express_in) {
  std::vector<std::string> expressions;
  expressions.emplace_back(ToMathematica(r3_element.x, express_in));
  expressions.emplace_back(ToMathematica(r3_element.y, express_in));
  expressions.emplace_back(ToMathematica(r3_element.z, express_in));
  return Apply("List", expressions);
}

template<typename D, typename... Qs>
std::string ToMathematica(Quantity<D> const& quantity,
                          ExpressIn<Qs...> express_in) {
  if constexpr (ExpressIn<Qs...>::is_default) {
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
                          ExpressIn<Qs...> express_in) {
  return ToMathematica(vector.coordinates(), express_in);
}

template<typename S, typename F, typename... Qs>
std::string ToMathematica(Bivector<S, F> const& bivector,
                          ExpressIn<Qs...> express_in) {
  return ToMathematica(bivector.coordinates(), express_in);
}

template<typename V, typename... Qs>
std::string ToMathematica(Point<V> const & point,
                          ExpressIn<Qs...> express_in) {
  return ToMathematica(point - Point<V>(), express_in);
}

template<typename F, typename... Qs>
std::string ToMathematica(DegreesOfFreedom<F> const& degrees_of_freedom,
                          ExpressIn<Qs...> express_in) {
  return Apply(
      "List",
      std::vector<std::string>{
          ToMathematica(degrees_of_freedom.position(), express_in),
          ToMathematica(degrees_of_freedom.velocity(), express_in)});
}

template<typename Tuple, typename, typename... Qs>
std::string ToMathematica(Tuple const& tuple, ExpressIn<Qs...> express_in) {
  std::vector<std::string> expressions;
  expressions.reserve(std::tuple_size_v<Tuple>);
  TupleHelper<std::tuple_size_v<Tuple>, Tuple, Qs...>::ToMathematicaStrings(
      tuple, expressions, express_in);
  return Apply("List", expressions);
}

template<typename R, typename, typename, typename... Qs>
std::string ToMathematica(R const ref,
                          ExpressIn<Qs...> express_in) {
  return Apply("List",
               std::vector<std::string>{
                   ToMathematica(ref.time, express_in),
                   ToMathematica(ref.degrees_of_freedom, express_in)});
}

template<typename... Qs>
std::string ToMathematica(
    astronomy::OrbitalElements::EquinoctialElements const& elements,
    ExpressIn<Qs...> express_in) {
  return ToMathematica(std::make_tuple((elements.t - J2000),
                                       elements.a,
                                       elements.h,
                                       elements.k,
                                       elements.λ,
                                       elements.p,
                                       elements.q,
                                       elements.pʹ,
                                       elements.qʹ),
                       express_in);
}

template<typename... Qs>
std::string ToMathematica(std::string const& str,
                          ExpressIn<Qs...>) {
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

}  // namespace internal_mathematica
}  // namespace mathematica
}  // namespace principia
