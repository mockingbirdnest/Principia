#pragma once

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
