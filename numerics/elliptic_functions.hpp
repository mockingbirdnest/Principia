#pragma once

#include "quantities/quantities.hpp"

// This code is a straightforward translation in C++ of:
// Fukushima, Toshio. (2012). xgscd.txt (Fortran program package to compute the
// Jacobian elliptic functions, sn(u|m), cn(u|m), dn(u|m)).
// Downloaded from:
// https://www.researchgate.net/publication/233903220_xgscdtxt_Fortran_program_package_to_compute_the_Jacobian_elliptic_functions_snum_cnum_dnum
namespace principia {
namespace numerics {
namespace internal_elliptic_functions {

using quantities::Angle;

Angle JacobiAmplitude(double u, double mc);

void JacobiSNCNDN(double u, double mc, double& s, double& c, double& d);

}  // namespace internal_elliptic_functions

using internal_elliptic_functions::JacobiAmplitude;
using internal_elliptic_functions::JacobiSNCNDN;

}  // namespace numerics
}  // namespace principia
