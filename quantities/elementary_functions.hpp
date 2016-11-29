
#pragma once

#include "quantities/quantities.hpp"

namespace principia {
namespace quantities {

// We need this for templates, for consistency with the dimensionful Sqrt.
// Equivalent to |std::sqrt(x)|.
double Sqrt(double x);
template<typename D>
SquareRoot<Quantity<D>> Sqrt(Quantity<D> const& x);

double Cbrt(double x);
template<typename D>
CubeRoot<Quantity<D>> Cbrt(Quantity<D> const& x);

double Sin(Angle const& α);
double Cos(Angle const& α);
double Tan(Angle const& α);

Angle ArcSin(double x);
Angle ArcCos(double x);
Angle ArcTan(double y, double x = 1);
template<typename D>
Angle ArcTan(Quantity<D> const& y, Quantity<D> const& x);

// We consider hyperbolic functions as dealing with quotients of arc length to
// curvature radius in the hyperbolic plane, which are angles. This explains the
// use of "arc" for inverse functions.

double Sinh(Angle const& α);
double Cosh(Angle const& α);
double Tanh(Angle const& α);

Angle ArcSinh(double x);
Angle ArcCosh(double x);
Angle ArcTanh(double x);
}  // namespace quantities
}  // namespace principia

#include "quantities/elementary_functions_body.hpp"
