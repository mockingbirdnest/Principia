#pragma once

#include "quantities/quantities.hpp"

// This code is derived from: [Fuk12a].  The original code has been translated
// into C++ and adapted to the needs of this project.
namespace principia {
namespace numerics {
namespace _elliptic_functions {
namespace internal {

using namespace principia::quantities::_quantities;

Angle JacobiAmplitude(Angle const& u, double mc);

void JacobiSNCNDN(Angle const& u, double mc, double& s, double& c, double& d);

}  // namespace internal

using internal::JacobiAmplitude;
using internal::JacobiSNCNDN;

}  // namespace _elliptic_functions
}  // namespace numerics
}  // namespace principia
