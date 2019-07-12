#pragma once

#include "quantities/quantities.hpp"

// This code is derived from: Fukushima, Toshio. (2018). xelbdj.txt: Fortran
// test driver for "elbdj"/"relbdj", subroutines to compute the double/single
// precision general incomplete elliptic integrals of all three kinds,
// DOI: 10.13140/RG.2.2.11113.80489, License: MIT.  Downloaded from:
// https://www.researchgate.net/publication/322702514_xelbdjtxt_Fortran_test_driver_for_elbdjrelbdj_subroutines_to_compute_the_doublesingle_precision_general_incomplete_elliptic_integrals_of_all_three_kinds
// The original code has been translated into C++ and adapted to the needs of
// this project.
namespace principia {
namespace numerics {
namespace internal_elliptic_integrals {

using quantities::Angle;

void FukushimaEllipticBDJ(Angle const& φ,
                          double n,
                          double mc,
                          Angle& b,
                          Angle& d,
                          Angle& j);

Angle EllipticE(Angle const& φ,
                double mc);

Angle EllipticF(Angle const& φ,
                double mc);

Angle EllipticΠ(Angle const& φ,
                double n,
                double mc);

void EllipticEFΠ(Angle const& φ,
                 double n,
                 double mc,
                 Angle& e,
                 Angle& f,
                 Angle& ᴨ);

Angle EllipticK(double mc);

}  // namespace internal_elliptic_integrals

using internal_elliptic_integrals::EllipticE;
using internal_elliptic_integrals::EllipticF;
using internal_elliptic_integrals::EllipticΠ;
using internal_elliptic_integrals::EllipticEFΠ;
using internal_elliptic_integrals::EllipticK;
using internal_elliptic_integrals::FukushimaEllipticBDJ;

}  // namespace numerics
}  // namespace principia
