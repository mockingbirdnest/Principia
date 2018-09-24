
#pragma once

#include "astronomy/fortran_astrodynamics_toolkit.hpp"

#include <algorithm>
#include <cmath>

namespace principia {
namespace astronomy {
namespace fortran_astrodynamics_toolkit {

using numerics::FixedVector;

template<int nmodel, int mmodel>
R3Element<double> ComputeGravityAccelerationLear(
    R3Element<double> const& rgr,
    double const mu,
    double const rbar,
    FixedMatrix<double, nmodel + 1, nmodel + 1> const& cnm,
    FixedMatrix<double, nmodel + 1, nmodel + 1> const& snm) {
  FixedMatrix<double, nmodel + 1, nmodel + 1> pnm;
  FixedMatrix<double, nmodel + 1, nmodel + 1> ppnm;
  FixedVector<double, nmodel + 1> cm;
  FixedVector<double, nmodel + 1> sm;
  FixedVector<double, nmodel + 1> pn;
  FixedVector<double, nmodel + 1> rb;
  FixedVector<double, nmodel + 1> ppn;
  R3Element<double> asph;
  double e1, e2, e3, e4, e5, r1, r2, t1, t3, absr, sphi, cphi, tcm, tsm, tsnm,
      tcnm, tpnm;
  int nm1, nm2;

  for (int n = 2; n <= nmodel; ++n) {
    pnm[n - 1][n] = 0.0;
  }

  e1 = rgr.x * rgr.x + rgr.y * rgr.y;
  r2 = e1 + rgr.z * rgr.z;
  absr = std::sqrt(r2);
  r1 = std::sqrt(e1);
  sphi = rgr.z / absr;
  cphi = r1 / absr;
  if (r1 == 0.0) {
    sm[1] = 0.0;
    cm[1] = 1.0;
  } else {
    sm[1] = rgr.y / r1;
    cm[1] = rgr.x / r1;
  }
  rb[1] = rbar / absr;
  rb[2] = rb[1] * rb[1];
  sm[2] = 2.0 * cm[1] * sm[1];
  cm[2] = 2.0 * cm[1] * cm[1] - 1.0;
  pn[1] = sphi;
  pn[2] = (3.0 * sphi * sphi - 1.0) / 2.0;
  ppn[1] = 1.0;
  ppn[2] = 3.0 * sphi;
  pnm[1][1] = 1.0;
  pnm[2][2] = 3.0 * cphi;
  pnm[2][1] = ppn[2];
  ppnm[1][1] = -sphi;
  ppnm[2][2] = -6.0 * sphi * cphi;
  ppnm[2][1] = 3.0 - 6.0 * sphi * sphi;
  if (nmodel >= 3) {
    for (int n = 3; n <= nmodel; ++n) {
      nm1 = n - 1;
      nm2 = n - 2;
      rb[n] = rb[nm1] * rb[1];
      sm[n] = 2.0 * cm[1] * sm[nm1] - sm[nm2];
      cm[n] = 2.0 * cm[1] * cm[nm1] - cm[nm2];
      e1 = 2 * n - 1;
      pn[n] = (e1 * sphi * pn[nm1] - nm1 * pn[nm2]) / n;
      ppn[n] = sphi * ppn[nm1] + n * pn[nm1];
      pnm[n][n] = e1 * cphi * pnm[nm1][nm1];
      ppnm[n][n] = -n * sphi * pnm[n][n];
    }
    for (int n = 3; n <= nmodel; ++n) {
      nm1 = n - 1;
      e1 = (2 * n - 1) * sphi;
      e2 = -n * sphi;
      for (int m = 1; m <= nm1; ++m) {
        e3 = pnm[nm1][m];
        e4 = n + m;
        e5 = (e1 * e3 - (e4 - 1.0) * pnm[n - 2][m]) / (n - m);
        pnm[n][m] = e5;
        ppnm[n][m] = e2 * e5 + e4 * e3;
      }
    }
  }

  asph.x = -1.0;
  asph.z = 0.0;

  for (int n = 2; n <= nmodel; ++n) {
    e1 = cnm[n][0] * rb[n];
    asph.x = asph.x - (n + 1) * e1 * pn[n];
    asph.z = asph.z + e1 * ppn[n];
  }
  asph.z = cphi * asph.z;
  t1 = 0.0;
  t3 = 0.0;
  asph.y = 0.0;

  for (int n = 2; n <= nmodel; ++n) {
    e1 = 0.0;
    e2 = 0.0;
    e3 = 0.0;
    for (int m = 1; m <= std::min(n, mmodel); ++m) {
      tsnm = snm[n][m];
      tcnm = cnm[n][m];
      tsm = sm[m];
      tcm = cm[m];
      tpnm = pnm[n][m];
      e4 = tsnm * tsm + tcnm * tcm;
      e1 = e1 + e4 * tpnm;
      e2 = e2 + m * (tsnm * tcm - tcnm * tsm) * tpnm;
      e3 = e3 + e4 * ppnm[n][m];
    }
    t1 = t1 + (n + 1) * rb[n] * e1;
    asph.y = asph.y + rb[n] * e2;
    t3 = t3 + rb[n] * e3;
  }

  e4 = mu / r2;
  asph.x = e4 * (asph.x - cphi * t1);
  asph.y = e4 * asph.y;
  asph.z = e4 * (asph.z + t3);

  e5 = asph.x * cphi - asph.z * sphi;

  R3Element<double> agr;
  agr.x = e5 * cm[1] - asph.y * sm[1];
  agr.y = e5 * sm[1] + asph.y * cm[1];
  agr.z = asph.x * sphi + asph.z * cphi;
  return agr;
}

}  // namespace fortran_astrodynamics_toolkit
}  // namespace astronomy
}  // namespace principia
