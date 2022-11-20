#pragma once

#include <array>

namespace principia {
namespace testing_utilities {

double GoldsteinPrice(double x₁, double x₂);

std::array<double, 2> GradGoldsteinPrice(double x₁, double x₂);

}  // namespace testing_utilities
}  // namespace principia
