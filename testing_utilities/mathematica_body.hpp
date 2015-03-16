#pragma once

#include "testing_utilities/mathematica.hpp"

#include <string>
#include <vector>

namespace principia {

using quantities::DebugString;

namespace testing_utilities {

template<typename T>
std::string MathematicaFunction(std::string const& function,
                                std::vector<T> const& arguments) {
  std::string result;
  result += function;
  result += "[";
  for (int i = 0; i != arguments.size(); ++i) {
    result += ToMathematica(arguments[i]);
    result += (i + 1 == arguments.size() ? "]" : ",");
  }
  return result;
}

inline std::string MathematicaFunction(
    std::string const& function,
    std::vector<std::string> const& arguments) {
  return MathematicaFunction<std::string>(function, arguments);
}

template<typename T>
std::string MathematicaOption(std::string const& name, T const& right) {
  return MathematicaFunction("Rule", {right});
}

template<typename T>
std::string MathematicaAssign(std::string const& name, T const& right) {
  return MathematicaFunction("Set", {name, ToMathematica(right)}) + ";\n";
}


template<typename T, typename U>
std::string MathematicaPlottableDataset(std::vector<T> x, std::vector<U> y) {
  std::vector<std::string> const xy = {ToMathematica(x), ToMathematica(y)};
  return MathematicaFunction("Transpose", {ToMathematica(xy)});
}

template<typename T>
std::string ToMathematica(std::vector<T> list) {
  return MathematicaFunction("List", list);
}

inline std::string ToMathematica(double const& real) {
  if (isinf(real)) {
    if (real > 0) {
      return "Infinity";
    } else {
      return MathematicaFunction("Minus", {"Infinity"});
    }
  } else if (isnan(real)) {
    return "Indeterminate";
  } else {
    std::string s = DebugString(real);
    s.replace(s.find("e"), 1, "*^");
    return MathematicaFunction("SetPrecision", {s, "MachinePrecision"});
  }
}

template<typename D>
std::string ToMathematica(Quantity<D> const& quantity) {
  std::string s = DebugString(quantity);
  s.replace(s.find("e"), 1, "*^");
  std::string number = ToMathematica(quantity / SIUnit<Quantity<D>>());
  std::size_t split = s.find(" ");
  std::string units = MathematicaEscape(s.substr(split, s.size()));
  return MathematicaFunction(
      "SetPrecision",
      {MathematicaFunction("Quantity", {number, units}),
       "MachinePrecision"});
}

// Wraps the string in quotes.
// TODO(egg): escape things properly.
inline std::string MathematicaEscape(std::string const& str) {
  std::string result = {"\""};
  result += str;
  result += "\"";
  return result;
}

inline std::string ToMathematica(std::string const& str) {
  return str;
}

}  // namespace testing_utilities
}  // namespace principia
