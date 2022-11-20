#include "testing_utilities/optimization_test_functions.hpp"

#include "quantities/elementary_functions.hpp"

namespace principia {
namespace testing_utilities {
namespace {
using quantities::Pow;
}  // namespace

double GoldsteinPrice(double const x₁, double const x₂) {
  return (1 + Pow<2>(x₁ + x₂ + 1) * (19 - 14 * x₁ + 3 * Pow<2>(x₁) - 14 * x₂ +
                                     6 * x₁ * x₂ + 3 * Pow<2>(x₂))) *
         (30 +
          Pow<2>(2 * x₁ - 3 * x₂) * (18 - 32 * x₁ + 12 * Pow<2>(x₁) + 48 * x₂ -
                                     36 * x₁ * x₂ + 27 * Pow<2>(x₂)));
}

std::array<double, 2> GradGoldsteinPrice(double const x₁, double const x₂) {
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

}  // namespace testing_utilities
}  // namespace principia
