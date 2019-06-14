#pragma once

// This code is a straightforward translation in C++ of:
// Fukushima, Toshio. (2012). xgscd.txt (Fortran program package to compute the
// Jacobian elliptic functions, sn(u|m), cn(u|m), dn(u|m)).
// Downloaded from:
// https://www.researchgate.net/publication/233903220_xgscdtxt_Fortran_program_package_to_compute_the_Jacobian_elliptic_functions_snum_cnum_dnum
namespace principia {
namespace numerics {

void Gscd(double const u, double const mc, double& s, double& c, double& d);

double Elk(double const mc);

}  // namespace numerics
}  // namespace principia
