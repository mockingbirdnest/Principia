#pragma once

#include <string>
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
// |x| and |y| should have the same size.
template<typename T, typename U>
quantities::Dimensionless PearsonProductMomentCorrelationCoefficient(
    std::vector<T> const& x,
    std::vector<U> const& y);

// The slope of the least-squares linear regression to the dataset
// with abscissae |x| and ordinates |y|. |x| and |y| should have the same size.
template<typename T, typename U>
quantities::Quotient<U, T> Slope(std::vector<T> const& x,
                                 std::vector<U> const& y);

// Mathematica input for a bidimensional dataset, copyable from the command line
// (includes a call to StringReplace to remove stray newlines). The resulting
// expression is a |List| of pairs {xᵢ, yᵢ}, which can be given as an argument
// to |ListPlot|. |x| and |y| should have the same size.
// The string contains one pair per line so as to make the command line output
// legible.
// Sample output:
// ToExpression[StringReplace["
// {{-6.9897000433601875e-001,-7.5318843669225339e+000},
//  {-7.4036268949424389e-001,-7.7377719257263520e+000},
//  {-7.8175537465246892e-001,-7.9437305531254143e+000},
//  {-8.2314805981069394e-001,-8.2105711052190316e+000},
//  {-8.6454074496891897e-001,-8.4686936817521428e+000},
//  {-9.0593343012714411e-001,-8.6385743611519761e+000},
//  {-9.4732611528536914e-001,-8.7871008196313465e+000},
//  {-9.8871880044359428e-001,-9.0271835254377528e+000},
//  {-1.0301114856018194e+000,-9.1953571177216453e+000},
//  {-1.0715041707600443e+000,-9.3888750926237385e+000},
//  {-1.1128968559182695e+000,-9.6011017004651631e+000},
//  {-1.1542895410764946e+000,-9.8238916688942766e+000},
//  {-1.1956822262347195e+000,-1.0045949826966634e+001},
//  {-1.2370749113929447e+000,-1.0252483352132446e+001},
//  {-1.2784675965511698e+000,-1.0429742060890906e+001},
//  {-1.3198602817093947e+000,-1.0640096857743393e+001},
//  {-1.3612529668676199e+000,-1.0850054861801775e+001},
//  {-1.4026456520258450e+000,-1.1063135945952322e+001},
//  {-1.4440383371840699e+000,-1.1247522304038885e+001},
//  {-1.4854310223422951e+000,-1.1453886857806401e+001},
//  {-1.5268237075005202e+000,-1.1666519487547756e+001},
//  {-1.5682163926587451e+000,-1.1885365812893951e+001}}",
// {"e"->"*^", "\n"->"", " "->""}]];
template<typename T, typename U>
std::string BidimensionalDatasetMathematicaInput(std::vector<T> const& x,
                                                 std::vector<U> const& y);

}  // namespace testing_utilities
}  // namespace principia

#include "testing_utilities/statistics_body.hpp"
