
#include "numerics/elliptic_functions.hpp"

#include "glog/logging.h"
#include "numerics/elliptic_integrals.hpp"
#include "quantities/numbers.hpp"

namespace principia {
namespace numerics {
namespace {

void Scd2(double u, double mc, double& s, double& c, double& d);

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
void Scd2(double const u, double const mc, double& s, double& c, double& d) {
  double m, uA, uT, u0, v, a, b, x, y, z, my, mc2, m2, xz, w;
  int n;  // TODO(phl): Used after the loop.

  constexpr double B10 = 1.0 / 24.0;
  constexpr double B11 = 1.0 / 6.0;
  constexpr double B20 = 1.0 / 720.0;
  constexpr double B21 = 11.0 / 180.0;
  constexpr double B22 = 1.0 / 45.0;

  m = 1.0 - mc;
  uA = 1.76269 + mc * 1.16357;
  uT = 5.217e-3 - m * 2.143e-3;
  u0 = u;

  for (n = 0; n <= 20; ++n) {
    if (u0 < uT) {
      break;
    }
    LOG_IF(FATAL, n == 20) << "(scd2) Too large input argument: u=" << u;
    u0 = u0 * 0.5;
  }
  v = u0 * u0;
  a = 1.0;
  b = v * (0.5 - v * (B10 + m * B11 - v * (B20 + m * (B21 + m * B22))));
  if (u < uA) {
    for (int j = 1; j <= n; ++j) {
      y = b * (a * 2.0 - b);
      z = a * a;
      my = m * y;
      b = (y * 2.0) * (z - my);
      a = z * z - my * y;
    }
  } else {
    for (int j = 1; j <= n; ++j) {
      y = b * (a * 2.0 - b);
      z = a * a;
      my = m * y;
      if (z < my * 2.0) {
        c = a - b;
        mc2 = mc * 2.0;
        m2 = m * 2.0;
        for (int i = j; i <= n; ++i) {
          x = c * c;
          z = a * a;
          w = m * x * x - mc * z * z;
          xz = x * z;
          c = mc2 * xz + w;
          a = m2 * xz - w;
        }
        c = c / a;
        x = c * c;
        s = sqrt(1.0 - x);
        d = sqrt(mc + m * x);
        return;
      }
      b = (y * 2.0) * (z - my);
      a = z * z - my * y;
    }
  }
  b = b / a;
  y = b * (2.0 - b);
  c = 1.0 - b;
  s = sqrt(y);
  d = sqrt(1.0 - m * y);
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
void Gscd(double const u, double const mc, double& s, double& c, double& d) {
  double m, kc, ux, k, kh, kh3, kh5, kh7, k2, k3, k4, sx;

  m = 1.0 - mc;
  kc = sqrt(mc);
  ux = abs(u);
  if (ux < 0.785) {
    Scd2(ux, mc, s, c, d);
  } else {
    k = Elk(mc);
    kh = k * 0.5;
    kh3 = k * 1.5;
    kh5 = k * 2.5;
    kh7 = k * 3.5;
    k2 = k * 2.0;
    k3 = k * 3.0;
    k4 = k * 4.0;
    ux = ux - k4 * static_cast<double>(static_cast<int>(ux / k4));
    if (ux < kh) {
      Scd2(ux, mc, s, c, d);
    } else if (ux < k) {
      ux = k - ux;
      Scd2(ux, mc, s, c, d);
      sx = c / d;
      c = kc * s / d;
      s = sx;
      d = kc / d;
    } else if (ux < kh3) {
      ux = ux - k;
      Scd2(ux, mc, s, c, d);
      sx = c / d;
      c = -kc * s / d;
      s = sx;
      d = kc / d;
    } else if (ux < k2) {
      ux = k2 - ux;
      Scd2(ux, mc, s, c, d);
      c = -c;
    } else if (ux < kh5) {
      ux = ux - k2;
      Scd2(ux, mc, s, c, d);
      s = -s;
      c = -c;
    } else if (ux < k3) {
      ux = k3 - ux;
      Scd2(ux, mc, s, c, d);
      sx = -c / d;
      c = -kc * s / d;
      s = sx;
      d = kc / d;
    } else if (ux < kh7) {
      ux = ux - k3;
      Scd2(ux, mc, s, c, d);
      sx = -c / d;
      c = kc * s / d;
      s = sx;
      d = kc / d;
    } else {
      ux = k4 - ux;
      Scd2(ux, mc, s, c, d);
      s = -s;
    }
  }
  if (u < 0.0) {
    s = -s;
  }
}

}  // namespace numerics
}  // namespace principia
