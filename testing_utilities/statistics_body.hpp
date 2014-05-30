#pragma once

#include <string>
#include <vector>

#include "quantities/dimensionless.hpp"
#include "quantities/elementary_functions.hpp"
#include "quantities/quantities.hpp"
#include "testing_utilities/numerics.hpp"

namespace principia {
namespace testing_utilities {

template<typename T>
T Mean(std::vector<T> const& x) {
  T total;
  for (T const& x_i : x) {
    total += x_i;
  }
  return total / x.size();
}

template<typename T>
quantities::Product<T, T> Variance(std::vector<T> const& x) {
  return Covariance(x, x);
}

template<typename T, typename U>
quantities::Product<T, U> Covariance(std::vector<T> const& x,
                                     std::vector<U> const& y) {
  std::vector<quantities::Product<T, U>> integrand(x.size());
  T mean_x = Mean(x);
  U mean_y = Mean(y);
  for (std::size_t i = 0; i < x.size(); ++i) {
    integrand[i] = (x[i] - mean_x) * (y[i] - mean_y);
  }
  return Mean(integrand);
}

template<typename T>
T StandardDeviation(std::vector<T> const& x) {
  return quantities::Sqrt(Variance(x));
}

template<typename T, typename U>
quantities::Dimensionless PearsonProductMomentCorrelationCoefficient(
    std::vector<T> const& x,
    std::vector<U> const& y) {
  return Covariance(x, y) / (StandardDeviation(x) * StandardDeviation(y));
}

template<typename T, typename U>
quantities::Quotient<U, T> Slope(std::vector<T> const& x,
                                 std::vector<U> const& y) {
  return Covariance(x, y) / Variance(x);
}

template<typename T, typename U>
std::string BidimensionalDatasetMathematicaInput(std::vector<T> const& x,
                                                 std::vector<U> const& y) {
  static std::string const mathematica_line =
      "(*****************************************************)";
  std::string result = mathematica_line + "\n";
  result += "ToExpression[StringReplace[\"\n{";
  for (std::size_t i = 0; i < x.size(); ++i) {
    result += "{";
    // We use |ToString(Dimensionless const&)| in order to get enough digits.
    result += quantities::ToString(DoubleValue(x[i]));
    result += ",";
    result += quantities::ToString(DoubleValue(y[i]));
    result += "}";
    if (i + 1 < x.size()) {
      result += ",\n ";
    }
  }
  result += "}\",\n{\"e\"->\"*^\", \"\\n\"->\"\", \" \"->\"\"}]];\n";
  result += mathematica_line;
  return result;
}

}  // namespace testing_utilities
}  // namespace principia
