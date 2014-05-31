#pragma once

#include<vector>

namespace principia {
namespace testing_utilities {

inline void ComputeHarmonicOscillatorForce(double const t,
                                              std::vector<double> const& q,
                                              std::vector<double>* result) {
  (*result)[0] = -q[0];
}

inline void ComputeHarmonicOscillatorVelocity(std::vector<double> const& p,
                                                 std::vector<double>* result) {
  (*result)[0] = p[0];
}

}  // namespace testing_utilities
}  // namespace principia
