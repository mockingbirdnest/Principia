
#include "numerics/elliptic_functions.hpp"

#include <tuple>

#include "glog/logging.h"
#include "numerics/combinatorics.hpp"
#include "numerics/elliptic_integrals.hpp"
#include "numerics/polynomial.hpp"
#include "numerics/polynomial_evaluators.hpp"
#include "quantities/elementary_functions.hpp"
#include "quantities/numbers.hpp"

namespace principia {

using quantities::Sqrt;

namespace numerics {
namespace {

void JacobiSNCNDNReduced(double u, double mc, double& s, double& c, double& d);

// Maclaurin series for Fukushima b₀.  These are polynomials in m that are used
// as coefficients of a polynomial in u₀².  The index gives the corresponding
// power of u₀².
PolynomialInMonomialBasis<double, double, 0, HornerEvaluator>
    fukushima_b₀_maclaurin_m_1(std::make_tuple(1.0 / 2.0));
PolynomialInMonomialBasis<double, double, 1, HornerEvaluator>
    fukushima_b₀_maclaurin_m_2(std::make_tuple(-1.0 / 24.0, -1.0 / 6.0));
PolynomialInMonomialBasis<double, double, 2, HornerEvaluator>
    fukushima_b₀_maclaurin_m_3(std::make_tuple(1.0 / 720.0,
                                               11.0 / 180.0,
                                               1.0 / 45.0));

// Double precision subroutine to compute three Jacobian elliptic functions
// simultaneously
//
//   For limited argument: 0 <= u < K/2
//
//     Reference: T. Fukushima, (2012) Numer. Math.
//       DOI 10.1007/s00211-012-0498-0
//       "Precise and Fast Computation of Jacobian Elliptic Functions by
//        Conditional Duplication"
//
//     Author: T. Fukushima Toshio.Fukushima@nao.ac.jp
//
//     Inputs: u = argument, mc = 1-m, 0 < mc <= 1
//
//     Output: s = sn(u|m), c=cn(u|m), d=dn(u|m)
//
void JacobiSNCNDNReduced(double const u,
                         double const mc,
                         double& s,
                         double& c,
                         double& d) {
  constexpr int max_reductions = 20;

  double const m = 1.0 - mc;
  double const uT = 5.217e-3 - 2.143e-3 * m;

  double u₀ = u;
  int n = 0;  // Note that this variable is used after the loop.
  for (; u₀ >= uT; ++n) {
    DCHECK_LE(n, max_reductions)
        << "u₀ = " << u₀ << " u = " << u << " mc = " << mc;
    u₀ = 0.5 * u₀;
  }

  double const b₀1 = fukushima_b₀_maclaurin_m_1.Evaluate(m);
  double const b₀2 = fukushima_b₀_maclaurin_m_2.Evaluate(m);
  double const b₀3 = fukushima_b₀_maclaurin_m_3.Evaluate(m);
  PolynomialInMonomialBasis<double, double, 3, HornerEvaluator>
      fukushima_b₀_maclaurin_u₀²_3(std::make_tuple(0.0, b₀1, b₀2, b₀3));
  double const u₀² = u₀ * u₀;

  // We use the subscript i to indicate variables that are computed as part of
  // the iteration (Fukushima uses subscripts n and N).  This avoids confusion
  // between c (the result) and cᵢ (the intermediate numerator of c).
  double bᵢ = fukushima_b₀_maclaurin_u₀²_3.Evaluate(u₀²);

  double const uA = 1.76269 + 1.16357 * mc;
  bool const may_have_cancellation = u > uA;
  double aᵢ = 1.0;
  for (int i = 0; i < n; ++i) {
    double const yᵢ = bᵢ * (2.0 * aᵢ - bᵢ);
    double const zᵢ = aᵢ * aᵢ;
    double const myᵢ = m * yᵢ;
    if (may_have_cancellation && zᵢ < 2.0 * myᵢ) {
      double cᵢ = aᵢ - bᵢ;
      double const two_mc = 2.0 * mc;
      double const two_m = 2.0 * m;
      for (; i < n; ++i) {
        double const xᵢ = cᵢ * cᵢ;
        double const zᵢ = aᵢ * aᵢ;
        double const wᵢ = m * xᵢ * xᵢ - mc * zᵢ * zᵢ;
        double const xᵢzᵢ = xᵢ * zᵢ;
        cᵢ = two_mc * xᵢzᵢ + wᵢ;
        aᵢ = two_m * xᵢzᵢ - wᵢ;
      }
      c = cᵢ / aᵢ;
      double const c² = c * c;
      s = Sqrt(1.0 - c²);
      d = Sqrt(mc + m * c²);
      return;
    }
    bᵢ = 2.0 * yᵢ * (zᵢ - myᵢ);
    aᵢ = zᵢ * zᵢ - myᵢ * yᵢ;
  }
  bᵢ = bᵢ / aᵢ;
  double const yᵢ = bᵢ * (2.0 - bᵢ);
  c = 1.0 - bᵢ;
  s = Sqrt(yᵢ);
  d = Sqrt(1.0 - m * yᵢ);
}

}  // namespace

// Double precision subroutine to compute three Jacobian elliptic functions
// simultaneously
//
//   For general argument: -infty < u < infty
//
//     Reference: T. Fukushima, (2012) Numer. Math.
//     DOI 10.1007/s00211-012-0498-0
//       "Precise and Fast Computation of Jacobian Elliptic Functions by
//        Conditional Duplication"
//
//     Author: T. Fukushima Toshio.Fukushima@nao.ac.jp
//
//     Used subprograms: scd2, elk
//
//     Inputs: u = argument, mc = 1-m, 0 < mc <= 1
//
//     Output: s = sn(u|m), c=cn(u|m), d=dn(u|m)
//
void JacobiSNCNDN(double const u,
                  double const mc,
                  double& s,
                  double& c,
                  double& d) {
  double m, kc, ux, k, kh, kh3, kh5, kh7, k2, k3, k4, sx;

  m = 1.0 - mc;
  kc = sqrt(mc);
  ux = abs(u);
  if (ux < 0.785) {
    JacobiSNCNDNReduced(ux, mc, s, c, d);
  } else {
    k = EllipticK(mc);
    kh = k * 0.5;
    kh3 = k * 1.5;
    kh5 = k * 2.5;
    kh7 = k * 3.5;
    k2 = k * 2.0;
    k3 = k * 3.0;
    k4 = k * 4.0;
    ux = ux - k4 * static_cast<double>(static_cast<int>(ux / k4));
    if (ux < kh) {
      JacobiSNCNDNReduced(ux, mc, s, c, d);
    } else if (ux < k) {
      ux = k - ux;
      JacobiSNCNDNReduced(ux, mc, s, c, d);
      sx = c / d;
      c = kc * s / d;
      s = sx;
      d = kc / d;
    } else if (ux < kh3) {
      ux = ux - k;
      JacobiSNCNDNReduced(ux, mc, s, c, d);
      sx = c / d;
      c = -kc * s / d;
      s = sx;
      d = kc / d;
    } else if (ux < k2) {
      ux = k2 - ux;
      JacobiSNCNDNReduced(ux, mc, s, c, d);
      c = -c;
    } else if (ux < kh5) {
      ux = ux - k2;
      JacobiSNCNDNReduced(ux, mc, s, c, d);
      s = -s;
      c = -c;
    } else if (ux < k3) {
      ux = k3 - ux;
      JacobiSNCNDNReduced(ux, mc, s, c, d);
      sx = -c / d;
      c = -kc * s / d;
      s = sx;
      d = kc / d;
    } else if (ux < kh7) {
      ux = ux - k3;
      JacobiSNCNDNReduced(ux, mc, s, c, d);
      sx = -c / d;
      c = kc * s / d;
      s = sx;
      d = kc / d;
    } else {
      ux = k4 - ux;
      JacobiSNCNDNReduced(ux, mc, s, c, d);
      s = -s;
    }
  }
  if (u < 0.0) {
    s = -s;
  }
}

}  // namespace numerics
}  // namespace principia
