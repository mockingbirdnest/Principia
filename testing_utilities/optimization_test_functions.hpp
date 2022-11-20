#pragma once

#include <array>

namespace principia {
namespace testing_utilities {

// See https://www.sfu.ca/~ssurjano/goldpr.html.
double GoldsteinPrice(double x₁, double x₂);
std::array<double, 2> GradGoldsteinPrice(double x₁, double x₂);

// See https://www.sfu.ca/~ssurjano/branin.html.
double Branin(double x₁, double x₂);
std::array<double, 2> GradBranin(double x₁, double x₂);

}  // namespace testing_utilities
}  // namespace principia
