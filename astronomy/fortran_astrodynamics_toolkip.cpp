#include "astronomy/fortran_astrodynamics_toolkit.hpp"

#include <cmath>

namespace principia {
namespace astronomy {
namespace fortran_astrodynamics_toolkit {

#define FINDEX(expression) ((expression) - 1)
#define FINDEX0(expression) (expression)

using numerics::FixedVector;

template<int nmodel, int mmodel>
R3Element<double> Grav(R3Element<double> const& rgr,
          double const mu,
          double const rbar,
          FixedMatrix<double, nmodel, nmodel + 1> const& cnm,
          FixedMatrix<double, nmodel, nmodel + 1> const& snm) {
  FixedMatrix<double, nmodel, nmodel> pnm;
  FixedMatrix<double, nmodel, nmodel> ppnm;
  FixedVector<double, nmodel> cm;
  FixedVector<double, nmodel> sm;
  FixedVector<double, nmodel> pn;
  FixedVector<double, nmodel> rb;
  FixedVector<double, nmodel> ppn;
  R3Element<double> asph;
  double e1,e2,e3,e4,e5,r1,r2,t1,t3,absr,sphi,cphi,tcm,tsm,tsnm,tcnm,tpnm;
  int n,nm1,nm2,m;

  for (int n = 2; n <= nmodel; ++n) {
    pnm[FINDEX(n - 1)][FINDEX(n)] = 0.0;
  }

  e1 = rgr.x * rgr.x + rgr.y * rgr.y;
  r2 = e1 + rgr.z * rgr.z;
  absr = std::sqrt(r2);
  r1 = std::sqrt(e1);
  sphi = rgr.z / absr;
  cphi = r1 / absr;
  if (r1 == 0.0) {
    sm[FINDEX(1)] = 0.0;
    cm[FINDEX(1)] = 1.0;
  } else {
    sm[FINDEX(1)] = rgr.y / r1;
    cm[FINDEX(1)] = rgr.x / r1;
  }
  rb[FINDEX(1)] = rbar / absr;
  rb[FINDEX(2)] = rb[FINDEX(1)] * rb[FINDEX(1)];
  sm[FINDEX(2)] = 2.0 * cm[FINDEX(1)] * sm[FINDEX(1)];
  cm[FINDEX(2)] = 2.0 * cm[FINDEX(1)] * cm[FINDEX(1)] - 1.0;
  pn[FINDEX(1)] = sphi;
  pn[FINDEX(2)] = (3.0 * sphi * sphi - 1.0) / 2.0;
  ppn[FINDEX(1)] = 1.0;
  ppn[FINDEX(2)] = 3.0 * sphi;
  pnm[FINDEX(1)][FINDEX(1)] = 1.0;
  pnm[FINDEX(2)][FINDEX(2)] = 3.0 * cphi;
  pnm[FINDEX(2)][FINDEX(1)] = ppn[FINDEX(2)];
  ppnm[FINDEX(1)][FINDEX(1)] = -sphi;
  ppnm[FINDEX(2)][FINDEX(2)] = -6.0 * sphi * cphi;
  ppnm[FINDEX(2)][FINDEX(1)] = 3.0 - 6.0 * sphi * sphi;
  if (nmodel >= 3) {
    for (int n = 3; n <= nmodel; ++n) {
      nm1 = n - 1;
      nm2 = n - 2;
      rb[FINDEX(n)] = rb[FINDEX(nm1)] * rb[FINDEX(1)];
      sm[FINDEX(n)] = 2.0 * cm[FINDEX(1)] * sm[FINDEX(nm1)] - sm[FINDEX(nm2)];
      cm[FINDEX(n)] = 2.0 * cm[FINDEX(1)] * cm[FINDEX(nm1)] - cm[FINDEX(nm2)];
      e1 = 2 * n - 1;
      pn[FINDEX(n)] = (e1 * sphi * pn[FINDEX(nm1)] - nm1 * pn[FINDEX(nm2)]) / n;
      ppn[FINDEX(n)] = sphi * ppn[FINDEX(nm1)] + n * pn[FINDEX(nm1)];
      pnm[FINDEX(n)][FINDEX(n)] = e1 * cphi * pnm[FINDEX(nm1)][FINDEX(nm1)];
      ppnm[FINDEX(n)][FINDEX(n)] = -n * sphi * pnm[FINDEX(n)][FINDEX(n)];
    }
    for (int n = 3; n <= model; ++n) {
      nm1 = n - 1;
      e1 = (2 * n - 1) * sphi;
      e2 = -n * sphi;
      for (int m = 1, m <= nm1; ++m) {
        e3 = pnm[FINDEX(nm1)][FINDEX(m)];
        e4 = n + m;
        e5 = (e1 * e3 - (e4 - 1.0) * pnm[FINDEX(n - 2)][FINDEX(m)]) / (n - m);
        pnm[FINDEX(n)][FINDEX(m)] = e5;
        ppnm[FINDEX(n)][FINDEX(m)] = e2 * e5 + e4 * e3;
      }
    }
  }

  asph.x = -1.0;
  asph.z = 0.0;

  for (int n = 2; n <= nmodel; ++n) {
    e1 = 0.0;
    e2 = 0.0;
    e3 = 0.0;
    for (int m = 1; m <= std::min(n, mmodel); ++m) {
      tsnm = snm[FINDEX(n)][FINDEX0(m)];
      tcnm = cnm[FINDEX(n)][FINDEX0(m)];
      tsm = sm[FINDEX(m)];
      tcm = cm[FINDEX(m)];
      tpnm = pnm[FINDEX(n)][FINDEX(m)];
      e4 = tsnm * tsm + tcnm * tcm;
      e1 = e1 + e4 * tpnm;
      e2 = e2 + m * (tsnm * tcm - tcnm * tsm) * tpnm;
      e3 = e3 + e4 * ppnm[FINDEX(n)][FINDEX(m)];
    }
    t1 = t1 + (n + 1) * rb[FINDEX(n)] * e1;
    asph.y = asph.y + rb[FINDEX(n)] * e2;
    t3 = t3 + rb[FINDEX(n)] * e3;
  }

  e4 = mu / r2;
  asph.x = e4 * (asph.x - cphi * t1);
  asph.y = e4 * asph.y;
  asph.z = e4 * (asph.z + t3);

  e5 = asph.x * cphi - asph.z * sphi;

  R3Element<double> agr;
  agr.x = e5 * cm[FINDEX(1)] - asph.y * sm[FINDEX(1)];
  agr.y = e5 * sm[FINDEX(1)] + asph.y * cm[FINDEX(1)];
  agr.z = asph.x * sphi + asph.z * cphi;
  return agr;
}

auto x = Grav<1, 1>(R3Element<double>(),
                    0.0,
                    0.0,
                    FixedMatrix<double, 1, 2>({1.0, 1.0}),
                    FixedMatrix<double, 1, 2>({1.0, 1.0}));

}  // namespace fortran_astrodynamics_toolkit
}  // namespace astronomy
}  // namespace principia
