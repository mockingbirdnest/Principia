#pragma once

#include <array>

namespace principia {
namespace testing_utilities {
namespace _optimization_test_functions {
namespace internal {

// See https://www.sfu.ca/~ssurjano/branin.html.
double Branin(double x₁, double x₂);
std::array<double, 2> 𝛁Branin(double x₁, double x₂);

// See https://www.sfu.ca/~ssurjano/goldpr.html.
double GoldsteinPrice(double x₁, double x₂);
std::array<double, 2> 𝛁GoldsteinPrice(double x₁, double x₂);

// See https://www.sfu.ca/~ssurjano/hart3.html.
double Hartmann3(double x₁, double x₂, double x₃);
std::array<double, 3> 𝛁Hartmann3(double x₁, double x₂, double x₃);

}  // namespace internal

using internal::𝛁Branin;
using internal::𝛁GoldsteinPrice;
using internal::𝛁Hartmann3;
using internal::Branin;
using internal::GoldsteinPrice;
using internal::Hartmann3;

}  // namespace _optimization_test_functions
}  // namespace testing_utilities
}  // namespace principia
