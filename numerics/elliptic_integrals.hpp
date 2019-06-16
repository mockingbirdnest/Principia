#pragma once

// This code is a straightforward translation in C++ of:
// Fukushima, Toshio. (2018). xelbdj.txt: Fortran test driver for
// "elbdj"/"relbdj", subroutines to compute the double/single precision general
// incomplete elliptic integrals of all three kinds,
// DOI: 10.13140/RG.2.2.11113.80489, License: MIT.
// Downloaded from:
// https://www.researchgate.net/publication/322702514_xelbdjtxt_Fortran_test_driver_for_elbdjrelbdj_subroutines_to_compute_the_doublesingle_precision_general_incomplete_elliptic_integrals_of_all_three_kinds
namespace principia {
namespace numerics {

void FukushimaEllipticBDJ(double φ,
                          double n,
                          double mc,
                          double& b,
                          double& d,
                          double& j);

double EllipticK(double mc);

}  // namespace numerics
}  // namespace principia
