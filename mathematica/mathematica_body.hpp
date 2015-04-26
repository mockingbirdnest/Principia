#pragma once

#include "mathematica/mathematica.hpp"

#include <cmath>
#include <string>
#include <vector>

namespace principia {

using quantities::DebugString;

namespace mathematica {

inline std::string Apply(
    std::string const& function,
    std::vector<std::string> const& arguments) {
  std::string result;
  result += function;
  result += "[";
  for (int i = 0; i < arguments.size(); ++i) {
    result += arguments[i];
    result += (i + 1 == arguments.size() ? "]" : ",");
  }
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
std::string PlottableDataset(std::vector<T> x, std::vector<U> y) {
  std::vector<std::string> const xy = {ToMathematica(x), ToMathematica(y)};
  return Apply("Transpose", {ToMathematica(xy)});
}

template<typename T>
std::string ToMathematica(std::vector<T> list) {
  std::vector<std::string> expressions;
  expressions.reserve(list.size());
  for (auto const& expression : list) {
    expressions.emplace_back(ToMathematica(expression));
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
    return Apply("SetPrecision", {s, "MachinePrecision"});
  }
}

template<typename D>
std::string ToMathematica(Quantity<D> const& quantity) {
  std::string s = DebugString(quantity);
  s.replace(s.find("e"), 1, "*^");
  std::string const number = ToMathematica(quantity / SIUnit<Quantity<D>>());
  std::size_t const split = s.find(" ");
  std::string const units = Escape(s.substr(split, s.size()));
  return Apply(
      "SetPrecision",
      {Apply("Quantity", {number, units}), "MachinePrecision"});
}

template<typename... Types>
std::string ToMathematica(std::tuple<Types...> const& tuple) {
  std::vector<std::string> expressions;
  expressions.reserve(sizeof...(Types));
  //TODO(phl): There has to be a better way...
  for (int i = 0; i < sizeof...(Types); ++i) {
    switch (i) {
    case 0:
      expressions.emplace_back(ToMathematica(std::get<0>(tuple)));
      break;
    case 1:
      expressions.emplace_back(ToMathematica(std::get<1>(tuple)));
      break;
    case 2:
      expressions.emplace_back(ToMathematica(std::get<2>(tuple)));
      break;
    case 3:
      expressions.emplace_back(ToMathematica(std::get<3>(tuple)));
      break;
    case 4:
      expressions.emplace_back(ToMathematica(std::get<4>(tuple)));
      break;
    }
  }
  return Apply("List", expressions);
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

}  // namespace mathematica
}  // namespace principia
