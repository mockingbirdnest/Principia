#include "testing_utilities/optimization_test_functions.hpp"

#include "numerics/fixed_arrays.hpp"
#include "quantities/elementary_functions.hpp"
#include "quantities/numbers.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace testing_utilities {

using numerics::FixedMatrix;
using numerics::FixedVector;
using quantities::Cos;
using quantities::Pow;
using quantities::Sin;
using quantities::si::Radian;

namespace branin_parameters {
  constexpr double a = 1;
  constexpr double b = 5.1 / (4 * Pow<2>(π));
  constexpr double c = 5 / π;
  constexpr double r = 6;
  constexpr double s = 10;
  constexpr double t = 1 / (8 * π);
}  // namespace branin_parameters

namespace hartmann_parameters {
constexpr FixedVector<double, 4> α({1.0, 1.2, 3.0, 3.2});
constexpr FixedMatrix<double, /*rows=*/4, /*columns=*/3> A({3.0, 10, 30,
                                                            0.1, 10, 35,
                                                            3.0, 10, 30,
                                                            0.1, 10, 35});
constexpr FixedMatrix<double, /*rows=*/4, /*columns=*/3> P(
    {3689e-4, 1170e-4, 2673e-4,
     4699e-4, 4387e-4, 7470e-4,
     1091e-4, 8732e-4, 5547e-4,
      381e-4, 5743e-4, 8828e-4});
}  // namespace hartmann_parameters

double Branin(double const x₁, double const x₂) {
  using namespace branin_parameters;
  return a * Pow<2>(x₂ - b * Pow<2>(x₁) + c * x₁ - r) +
         s * (1 - t) * Cos(x₁ * Radian) + s;
}

std::array<double, 2> 𝛁Branin(double const x₁, double const x₂) {
  using namespace branin_parameters;
  double const g₁ = 2 * a * (c - 2 * b * x₁) * (-r + x₁ * (c - b * x₁) + x₂) +
                    s * (-1 + t) * Sin(x₁ * Radian);
  double const g₂ = 2 * a * (-r + x₁ * (c - b * x₁) + x₂);
  return {g₁, g₂};
}

double GoldsteinPrice(double const x₁, double const x₂) {
  return (1 + Pow<2>(x₁ + x₂ + 1) * (19 - 14 * x₁ + 3 * Pow<2>(x₁) - 14 * x₂ +
                                     6 * x₁ * x₂ + 3 * Pow<2>(x₂))) *
         (30 +
          Pow<2>(2 * x₁ - 3 * x₂) * (18 - 32 * x₁ + 12 * Pow<2>(x₁) + 48 * x₂ -
                                     36 * x₁ * x₂ + 27 * Pow<2>(x₂)));
}

std::array<double, 2> 𝛁GoldsteinPrice(double const x₁, double const x₂) {
  double const g₁ =
      24 * (-1 + 2 * x₁ - 3 * x₂) * (2 * x₁ - 3 * x₂) *
          (2 * x₁ - 3 * (1 + x₂)) *
          (1 +
           Pow<2>(1 + x₁ + x₂) * (19 + 3 * Pow<2>(x₁) + x₂ * (-14 + 3 * x₂) +
                                  2 * x₁ * (-7 + 3 * x₂))) +
      12 * (-2 + x₁ + x₂) * (-1 + x₁ + x₂) * (1 + x₁ + x₂) *
          (30 + Pow<2>(2 * x₁ - 3 * x₂) *
                    (18 + 12 * Pow<2>(x₁) - 4 * x₁ * (8 + 9 * x₂) +
                     3 * x₂ * (16 + 9 * x₂)));
  double const g₂ =
      -36 * (-1 + 2 * x₁ - 3 * x₂) * (2 * x₁ - 3 * x₂) *
          (2 * x₁ - 3 * (1 + x₂)) *
          (1 +
           Pow<2>(1 + x₁ + x₂) * (19 + 3 * Pow<2>(x₁) + x₂ * (-14 + 3 * x₂) +
                                  2 * x₁ * (-7 + 3 * x₂))) +
      12 * (-2 + x₁ + x₂) * (-1 + x₁ + x₂) * (1 + x₁ + x₂) *
          (30 + Pow<2>(2 * x₁ - 3 * x₂) *
                    (18 + 12 * Pow<2>(x₁) - 4 * x₁ * (8 + 9 * x₂) +
                     3 * x₂ * (16 + 9 * x₂)));
  return {g₁, g₂};
}

double Hartmann3(double const x₁, double const x₂, double const x₃) {
  using namespace hartmann_parameters;
  std::array<double, 3> const x{x₁, x₂, x₃};
  double result = 0;
  for (int i = 0; i < 4; ++i) {
    double exponent = 0;
    for (int j = 0; j < 3; ++j) {
      exponent -= A(i, j) * Pow<2>(x[j] - P(i, j));
    }
    result -= α[i] * std::exp(exponent);
  }
  return result;
}

std::array<double, 3> 𝛁Hartmann3(double const x₁,
                                 double const x₂,
                                 double const x₃) {
  using namespace hartmann_parameters;
  std::array<double, 3> const x{x₁, x₂, x₃};
  auto component =
      [&x](std::int64_t const k) {
        double result = 0;
        for (int i = 0; i < 4; ++i) {
          double exponent = 0;
          for (int j = 0; j < 3; ++j) {
            exponent -= A(i, j) * Pow<2>(x[j] - P(i, j));
          }
          result += 2 * α[i] * std::exp(exponent) * A(i, k) * (x[k] - P(i, k));
        }
        return result;
      };

  return {component(0), component(1), component(2)};
}


}  // namespace testing_utilities
}  // namespace principia
