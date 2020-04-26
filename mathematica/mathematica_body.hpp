
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
template<int index, typename Tuple, typename OptionalExpressIn>
struct TupleHelper : not_constructible {
  static void ToMathematicaStrings(Tuple const& tuple,
                                   std::vector<std::string>& expressions,
                                   OptionalExpressIn express_in) {
    TupleHelper<index - 1, Tuple, OptionalExpressIn>::ToMathematicaStrings(
        tuple, expressions, express_in);
    expressions.push_back(ToMathematica(std::get<index - 1>(tuple),
                          express_in));
  }
};

template<typename Tuple, typename OptionalExpressIn>
struct TupleHelper<0, Tuple, OptionalExpressIn> : not_constructible {
  static void ToMathematicaStrings(Tuple const& tuple,
                                   std::vector<std::string>& expressions,
                                   OptionalExpressIn express_in) {}
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

template<typename T, typename OptionalExpressIn>
std::string Option(std::string const& name,
                   T const& right,
                   OptionalExpressIn express_in) {
  return Apply("Rule", {name, ToMathematica(right, express_in)});
}

template<typename T, typename OptionalExpressIn>
std::string Assign(std::string const& name,
                   T const& right,
                   OptionalExpressIn express_in) {
  return Apply("Set", {name, ToMathematica(right, express_in)}) + ";\n";
}

template<typename T, typename U, typename OptionalExpressIn>
std::string PlottableDataset(std::vector<T> const& x,
                             std::vector<U> const& y,
                             OptionalExpressIn express_in) {
  std::vector<std::string> const xy = {ToMathematica(x, express_in),
                                       ToMathematica(y, express_in)};
  return Apply("Transpose", {ToMathematica(xy, express_in)});
}

template<typename T, typename OptionalExpressIn>
std::string ToMathematica(std::vector<T> const& list,
                          OptionalExpressIn express_in) {
  std::vector<std::string> expressions;
  expressions.reserve(list.size());
  for (auto const& expression : list) {
    expressions.emplace_back(ToMathematica(expression, express_in));
  }
  return Apply("List", expressions);
}

template<typename It, typename OptionalExpressIn>
std::string ToMathematica(It const begin, It const end,
                          OptionalExpressIn express_in) {
  std::vector<std::string> expressions;
  for (auto it = begin; it != end; ++it) {
    expressions.emplace_back(ToMathematica(*it, express_in));
  }
  return Apply("List", expressions);
}

template<typename OptionalExpressIn>
std::string ToMathematica(double const real,
                          OptionalExpressIn /*express_in*/) {
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

template<typename OptionalExpressIn>
std::string ToMathematica(Quaternion const& quaternion,
                          OptionalExpressIn /*express_in*/) {
  return Apply("Quaternion",
               {ToMathematica(quaternion.real_part()),
                ToMathematica(quaternion.imaginary_part().x),
                ToMathematica(quaternion.imaginary_part().y),
                ToMathematica(quaternion.imaginary_part().z)});
}

template<typename T, int size, typename OptionalExpressIn>
std::string ToMathematica(FixedVector<T, size> const& fixed_vector,
                          OptionalExpressIn express_in) {
  std::vector<std::string> expressions;
  for (int i = 0; i < size; ++i) {
    expressions.emplace_back(ToMathematica(fixed_vector[i], express_in));
  }
  return Apply("List", expressions);
}

template<typename T, typename OptionalExpressIn>
std::string ToMathematica(R3Element<T> const& r3_element,
                          OptionalExpressIn express_in) {
  std::vector<std::string> expressions;
  expressions.emplace_back(ToMathematica(r3_element.x, express_in));
  expressions.emplace_back(ToMathematica(r3_element.y, express_in));
  expressions.emplace_back(ToMathematica(r3_element.z, express_in));
  return Apply("List", expressions);
}

template<typename D, typename... Qs>
std::string ToMathematica(Quantity<D> const& quantity,
                          ExpressIn<Qs...> express_in) {
  return ToMathematica(express_in(quantity));
}

template<typename D, typename OptionalExpressIn>
std::string ToMathematica(Quantity<D> const& quantity,
                          std::nullopt_t express_in) {
  std::string s = DebugString(quantity);
  std::string const number = ToMathematica(quantity / si::Unit<Quantity<D>>);
  std::size_t const split = s.find(" ");
  std::string const units = Escape(s.substr(split, s.size()));
  return Apply("Quantity", {number, units});
}

template<typename S, typename F, typename OptionalExpressIn>
std::string ToMathematica(Vector<S, F> const& vector,
                          OptionalExpressIn express_in) {
  return ToMathematica(vector.coordinates(), express_in);
}

template<typename S, typename F, typename OptionalExpressIn>
std::string ToMathematica(Bivector<S, F> const& bivector,
                          OptionalExpressIn express_in) {
  return ToMathematica(bivector.coordinates(), express_in);
}

template<typename V, typename OptionalExpressIn>
std::string ToMathematica(Point<V> const & point,
                          OptionalExpressIn express_in) {
  return ToMathematica(point - Point<V>(), express_in);
}

template<typename F, typename OptionalExpressIn>
std::string ToMathematica(DegreesOfFreedom<F> const& degrees_of_freedom,
                          OptionalExpressIn express_in) {
  return Apply(
      "List",
      std::vector<std::string>{
          ToMathematica(degrees_of_freedom.position(), express_in),
          ToMathematica(degrees_of_freedom.velocity(), express_in)});
}

template<typename Tuple, typename, typename OptionalExpressIn>
std::string ToMathematica(Tuple const& tuple, OptionalExpressIn express_in) {
  std::vector<std::string> expressions;
  expressions.reserve(std::tuple_size_v<Tuple>);
  TupleHelper<std::tuple_size_v<Tuple>, Tuple, OptionalExpressIn>::
      ToMathematicaStrings(tuple, expressions, express_in);
  return Apply("List", expressions);
}

template<typename R, typename, typename, typename OptionalExpressIn>
std::string ToMathematica(R const ref,
                          OptionalExpressIn express_in) {
  return Apply("List",
               std::vector<std::string>{
                   ToMathematica(ref.time, express_in),
                   ToMathematica(ref.degrees_of_freedom, express_in)});
}

template<typename OptionalExpressIn>
std::string ToMathematica(
    astronomy::OrbitalElements::EquinoctialElements const& elements,
    OptionalExpressIn express_in) {
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

template<typename OptionalExpressIn>
std::string ToMathematica(std::string const& str,
                          OptionalExpressIn /*express_in*/) {
  return str;
}

inline std::string Escape(std::string const& str) {
  std::string result = "\"";
  for (const char c : str) {
    switch (c) {
      case '"':
        result += "\\\"";
        break;
      case '\\':
        result += "\\\\";
        break;
      default:
        result += c;
        break;
    }
  }
  result += "\"";
  return result;
}

inline Logger::Logger(std::filesystem::path const& path) : file_(path) {}

inline Logger::~Logger() {
  for (auto const& [name, values] : names_and_values_) {
    file_ << Assign(name, values);
  }
}

template<typename... Args>
void Logger::Append(std::string const& name, Args... args) {
  names_and_values_[name].push_back(ToMathematica(args...));
}

}  // namespace internal_mathematica
}  // namespace mathematica
}  // namespace principia
