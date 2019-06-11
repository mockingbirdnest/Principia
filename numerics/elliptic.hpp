#pragma once

// This code is a straightforward translation in C++ of:
// Fukushima, Toshio. (2018). xelbdj.txt: Fortran test driver for
// "elbdj"/"relbdj", subroutines to compute the double/single precision general
// incomplete elliptic integrals of all three kinds.
// DOI: 10.13140/RG.2.2.11113.80489.
// Downloaded from:
// https://www.researchgate.net/publication/322702514_xelbdjtxt_Fortran_test_driver_for_elbdjrelbdj_subroutines_to_compute_the_doublesingle_precision_general_incomplete_elliptic_integrals_of_all_three_kinds
namespace principia {
namespace numerics {

void Elbdj(double const phi,
           double const phic,
           double const n,
           double const mc,
           double& b,
           double& d,
           double& j);

double Cel(double const kc0,
           double const nc,
           double const aa,
           double const bb,
           int& err) ;

void Celbd(double const mc,double& elb, double& eld);

void Celbdj(double const nc,
            double const mc,
            double& bc,
            double& dc,
            double& jc);

void Elcbdj(double const c0,
            double const n,
            double const mc,
            double& b,
            double& d,
            double& j);

void Elsbdj(double const s0,
            double const n,
            double const mc,
            double& b,
            double& d,
            double& j);

void Serbd(double const y, double const m, double& b, double& d);

double Serj(double const y, double const n, double const m);

double Uatan(double const t, double const h);

}  // namespace numerics
}  // namespace principia
