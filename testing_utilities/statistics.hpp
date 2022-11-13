#pragma once

#include <string>
#include <vector>

#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"

// Various statistics on finite populations stored as |std::vector|s of
// |Quantity| or |Dimensionless|.

namespace principia {
namespace testing_utilities {
namespace internal_statistics {

using quantities::Product;
using quantities::Quotient;

// The population mean μ(x) = E[x].
template<typename T>
T Mean(std::vector<T> const& x);

// The population variance σ(x)² = Var[x] = E[x² - E[x]²].
// Computed as |Covariance(x, x)|.
template<typename T>
Product<T, T> Variance(std::vector<T> const& x);

// The population covariance Cov[x, y] = E[(x - E[x])(y - E[y])].
template<typename T, typename U>
Product<T, U> Covariance(std::vector<T> const& x, std::vector<U> const& y);

// The population standard deviation σ(x) = √(Var[x]).
template<typename T>
T StandardDeviation(std::vector<T> const& x);

// The Pearson product-moment correlation coefficient
// ρ(x, y) = Cov[x, y] / (σ(x) σ(y)).
// |x| and |y| should have the same size.
template<typename T, typename U>
double PearsonProductMomentCorrelationCoefficient(std::vector<T> const& x,
                                                  std::vector<U> const& y);

// The slope of the least-squares linear regression to the dataset
// with abscissae |x| and ordinates |y|. |x| and |y| should have the same size.
template<typename T, typename U>
Quotient<U, T> Slope(std::vector<T> const& x, std::vector<U> const& y);

// Mathematica input for a bidimensional dataset, copyable from the command line
// (includes a call to StringReplace to remove stray newlines). The resulting
// expression is a |List| of pairs {xᵢ, yᵢ}, which can be given as an argument
// to |ListPlot|. |x| and |y| should have the same size.
// The result contains one pair per line and is delimited by Mathematica
// comments so as to make the command line output more legible. It is not
// terminated by a newline.
// Sample output:
// (*****************************************************)
// ToExpression[StringReplace["
// {{-6.9897000433601875e-001,-7.5318843669225339e+000},
//  {-1.4854310223422951e+000,-1.1453886857806401e+001},
//  {-1.5268237075005202e+000,-1.1666519487547756e+001},
//  {-1.5682163926587451e+000,-1.1885365812893951e+001}}",
// {"e"->"*^", "\n"->"", " "->""}]];
// (*****************************************************)
template<typename T, typename U>
std::string BidimensionalDatasetMathematicaInput(std::vector<T> const& x,
                                                 std::vector<U> const& y);

}  // namespace internal_statistics

using internal_statistics::BidimensionalDatasetMathematicaInput;
using internal_statistics::Covariance;
using internal_statistics::Mean;
using internal_statistics::PearsonProductMomentCorrelationCoefficient;
using internal_statistics::Slope;
using internal_statistics::StandardDeviation;
using internal_statistics::Variance;

}  // namespace testing_utilities
}  // namespace principia

#include "testing_utilities/statistics_body.hpp"
