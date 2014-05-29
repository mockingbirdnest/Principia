#pragma once

#include <vector>

#include "quantities/dimensionless.hpp"
#include "quantities/quantities.hpp"

// Various statistics on finite populations stored as |std::vector|s of
// |Quantity| or |Dimensionless|.

namespace principia {
namespace testing_utilities {

// The population mean μ(x) = E[x].
template<typename T>
T Mean(std::vector<T> const& x);

// The population variance σ(x)² = Var[x] = E[x² - E[x]²].
// Computed as |Covariance(x, x)|.
template<typename T>
quantities::Product<T, T> Variance(std::vector<T> const& x);

// The population covariance Cov[x, y] = E[(x - E[x])(y - E[y])].
template<typename T, typename U>
quantities::Product<T, U> Covariance(std::vector<T> const& x,
                                     std::vector<U> const& y);

// The population standard deviation σ(x) = √(Var[x]).
template<typename T>
T StandardDeviation(std::vector<T> const& x);

// The Pearson product-moment correlation coefficient
// ρ(x, y) = Cov[x, y] / (σ(x) σ(y)).
// |x| and |y| should have the same |size()|.
template<typename T, typename U>
quantities::Dimensionless PearsonProductMomentCorrelationCoefficient(
    std::vector<T> const& x,
    std::vector<U> const& y);

// The slope of the least-squares linear regression to the dataset
// with abscissae |x| and ordinates |y|. |x| and |y| should have the same
// |size()|.
template<typename T, typename U>
quantities::Quotient<U, T> Slope(std::vector<T> const& x,
                                 std::vector<U> const& y);

}  // testing_utilities
}  // principia

#include "testing_utilities/statistics_body.hpp"
