#pragma once

#include "quantities/quantities.hpp"

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

// Same as above, but does not compute J.
void FukushimaEllipticBD(Angle const& φ, double mc, Angle& B_φǀm, Angle& D_φǀm);

// Returns the incomplete elliptic integral of the first kind F(φ|m), where
// m = 1 - mc.
Angle EllipticF(Angle const& φ, double mc);

// Returns the incomplete elliptic integral of the second kind E(φ|m), where
// m = 1 - mc.
Angle EllipticE(Angle const& φ, double mc);

// Returns the incomplete elliptic integral of the third kind Π(φ, n|m), where
// m = 1 - mc.
Angle EllipticΠ(Angle const& φ, double n, double mc);

// Computes the incomplete elliptic integrals the first and second kinds F(φ|m)
// and E(φ|m), where m = 1 - mc.
void EllipticFE(Angle const& φ, double mc, Angle& F_φǀm, Angle& E_φǀm);

// Computes the incomplete elliptic integrals of all three kinds F(φ|m), E(φ|m),
// and Π(φ, n|m), where m = 1 - mc.
void EllipticFEΠ(Angle const& φ,
                 double n,
                 double mc,
                 Angle& F_φǀm,
                 Angle& E_φǀm,
                 Angle& Π_φ_nǀm);

// Returns the complete elliptic integral of the first kind K(m), where
// m = 1 - mc.
Angle EllipticK(double mc);

}  // namespace internal_elliptic_integrals

using internal_elliptic_integrals::EllipticE;
using internal_elliptic_integrals::EllipticF;
using internal_elliptic_integrals::EllipticFE;
using internal_elliptic_integrals::EllipticFEΠ;
using internal_elliptic_integrals::EllipticK;
using internal_elliptic_integrals::EllipticΠ;
using internal_elliptic_integrals::FukushimaEllipticBD;
using internal_elliptic_integrals::FukushimaEllipticBDJ;

}  // namespace numerics
}  // namespace principia
