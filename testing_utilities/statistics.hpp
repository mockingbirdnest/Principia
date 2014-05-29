#pragma once

#include <vector>

#include "quantities/quantities.hpp"

namespace principia {
namespace testing_utilities {

template<typename T>
T Mean(std::vector<T> x);

template<typename T>
quantities::Product<T, T> Variance(std::vector<T> x);

template<typename T, typename U>
quantities::Product<T, U> Covariance(std::vector<T> x, std::vector<U> y);

template<typename T>
T StandardDeviation(std::vector<T> x);

template<typename T, typename U>
Dimensionless PearsonProductMomentCorrelation(std::vector<T> x,
                                              std::vector<U> y);

}  // testing_utilities
}  // principia

#include "testing_utilities/statistics_body.hpp"
