#pragma once

#include "quantities/quantities.hpp"

// Bibliography:
// [JE33] Jahnke and Emde (1933), Funktionentafeln mit Formeln und Kurven—Tables
// of functions with formulæ and curves.
// [JE38]  Jahnke and Emde (1938), Funktionentafeln mit Formeln und Kurven—Tables
// of functions with formulæ and curves.
// [JEL60] Jahnke, Emde, and Lösch (1960) Tafeln Höherer Funktionen—Tables of
// higher functions.
// [Bul65] Bulirsch (1965), Numerical Calculation of Elliptic Integrals and
// Elliptic Fuctions.
// [Bul69] Bulirsch (1969), Numerical Calculation of Elliptic Integrals and
// Elliptic Fuctions.  III.
// [Fuk11a] Fukushima (2011), Precise and fast computation of the general
// complete elliptic integral of the second kind.
// [Fuk11b] Fukushima (2011), Precise and fast computation of a general
// incomplete elliptic integral of second kind by half and double argument
// transformations.
// [Fuk12] Fukushima (2012), Precise and fast computation of a general
// incomplete elliptic integral of third kind by half and double argument
// transformations.
// [Fuk18] Fukushima (2018), xelbdj.txt: Fortran test driver for
// “elbdj”/“relbdj”, subroutines to compute the double / single precision
// general incomplete elliptic integrals of all three kinds.
// [OLBC10] Olver, Lozier, Boisvert, Clark Eds. (2010), NIST Handbook of
// Mathematical Functions.

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
// — the incomplete B is introduced by [Bul65], generalizing [JEL60];
// — the name “associate elliptic integral of the 2nd kind” is from [Fuk11b];
// — the integral J is introduced in [Fuk11c].
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

// Computes the incomplete elliptic integrals of all three kinds F(φ|m), E(φ|m),
// and Π(φ, n|m), where m = 1 - mc.
void EllipticFEΠ(Angle const& φ,
                 double n,
                 double mc,
                 Angle& F_φǀm,
                 Angle& E_φǀm,
                 Angle& Π_φ_nǀm);

Angle EllipticK(double mc);

}  // namespace internal_elliptic_integrals

using internal_elliptic_integrals::EllipticE;
using internal_elliptic_integrals::EllipticF;
using internal_elliptic_integrals::EllipticΠ;
using internal_elliptic_integrals::EllipticFEΠ;
using internal_elliptic_integrals::EllipticK;
using internal_elliptic_integrals::FukushimaEllipticBDJ;

}  // namespace numerics
}  // namespace principia
