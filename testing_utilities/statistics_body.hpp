#pragma once

#include <vector>

#include "quantities/elementary_functions.hpp"
#include "quantities/quantities.hpp"

namespace principia {
namespace testing_utilities {

template<typename T>
T Mean(std::vector<T> const& x) {
  T total;
  for(T const& x_i : x) {
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
  for(std::size_t i = 0; i < x.size(); ++i) {
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

}  // testing_utilities
}  // principia
