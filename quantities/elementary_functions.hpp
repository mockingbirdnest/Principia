
#pragma once

#include "quantities/quantities.hpp"

namespace principia {
namespace quantities {
// elementary_functions shares its internal namespace with quantities, because
// friendships otherwise become impossible to untangle.
namespace internal_quantities {

// Equivalent to |std::abs(x)|.
double Abs(double x);
template<typename D>
Quantity<D> Abs(Quantity<D> const& x);

// We need this for templates, for consistency with the dimensionful Sqrt.
// Equivalent to |std::sqrt(x)|.
double Sqrt(double x);
template<typename D>
SquareRoot<Quantity<D>> Sqrt(Quantity<D> const& x);

double Cbrt(double x);
template<typename D>
CubeRoot<Quantity<D>> Cbrt(Quantity<D> const& x);

// Equivalent to |std::pow(x, exponent)| unless -3 ≤ x ≤ 3, in which case
// explicit specialization yields multiplications statically.
template<int exponent>
constexpr double Pow(double x);
template<int exponent, typename D>
constexpr Exponentiation<Quantity<D>, exponent> Pow(Quantity<D> const& x);

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

}  // namespace internal_quantities

using internal_quantities::Abs;
using internal_quantities::ArcCos;
using internal_quantities::ArcCosh;
using internal_quantities::ArcSin;
using internal_quantities::ArcSinh;
using internal_quantities::ArcTan;
using internal_quantities::ArcTanh;
using internal_quantities::Cbrt;
using internal_quantities::Cos;
using internal_quantities::Cosh;
using internal_quantities::Pow;
using internal_quantities::Sin;
using internal_quantities::Sinh;
using internal_quantities::Sqrt;
using internal_quantities::Tan;
using internal_quantities::Tanh;

}  // namespace quantities
}  // namespace principia

#include "quantities/elementary_functions_body.hpp"
