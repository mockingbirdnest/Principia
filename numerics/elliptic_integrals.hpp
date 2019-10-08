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

// Computes the associate incomplete elliptic integrals of the second kind
// B(φ|m) and D(φ|m), as well as Fukushima’s associate incomplete elliptic
// integral of the third kind J(φ, n|m), where m = 1 - mc.
// NOTE(egg): As far as I can tell, the origins of the notation are as follows:
// — the integral D (complete and incomplete) is introduced in [JE33];
// — the complete integral B is introduced in the re-edition [JE38];
// — the incomplete B is introduced by [Buli65], generalizing [JEL60];
// — the name “associate elliptic integral of the 2nd kind” is from [Fuku11b];
// — the integral J is introduced in [Fuku11c].
void FukushimaEllipticBDJ(Angle const& φ,
                          double n,
                          double mc,
                          Angle& B_φǀm,
                          Angle& D_φǀm,
                          Angle& J_φ_nǀm);

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
