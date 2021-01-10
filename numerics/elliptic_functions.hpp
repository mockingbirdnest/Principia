#pragma once

#include "quantities/quantities.hpp"

// This code is derived from: [Fuk12a].  The original code has been translated
// into C++ and adapted to the needs of this project.
namespace principia {
namespace numerics {
namespace internal_elliptic_functions {

using quantities::Angle;

Angle JacobiAmplitude(Angle const& u, double mc);

void JacobiSNCNDN(Angle const& u, double mc, double& s, double& c, double& d);

}  // namespace internal_elliptic_functions

using internal_elliptic_functions::JacobiAmplitude;
using internal_elliptic_functions::JacobiSNCNDN;

}  // namespace numerics
}  // namespace principia
