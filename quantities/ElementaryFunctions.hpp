#pragma once

#include "quantities/Dimensionless.hpp"
#include "quantities/Quantities.hpp"

namespace principia {
namespace quantities {
Dimensionless Sqrt(Dimensionless const& x);
template<typename D>
SquareRoot<Quantity<D>> Sqrt(Quantity<D> const& x);
Dimensionless Log(Dimensionless const& x);
Dimensionless Log2(Dimensionless const& x);
Dimensionless Log10(Dimensionless const& x);
Dimensionless Exp(Dimensionless const& x);

Dimensionless Sin(Angle const& α);
Dimensionless Cos(Angle const& α);
Dimensionless Tan(Angle const& α);

Angle ArcSin(Dimensionless const& x);
Angle ArcCos(Dimensionless const& x);
Angle ArcTan(Dimensionless const& y, Dimensionless const& x = 1);
template<typename D>
Angle ArcTan(Quantity<D> const& y, Quantity<D> const& x);

// We consider hyperbolic functions as dealing with quotients of arc length to
// curvature radius in the hyperbolic plane, which are angles. This explains the
// use of "arc" for inverse functions.

Dimensionless Sinh(Angle const& α);
Dimensionless Cosh(Angle const& α);
Dimensionless Tanh(Angle const& α);

Angle ArcSinh(Dimensionless const& x);
Angle ArcCosh(Dimensionless const& x);
Angle ArcTanh(Dimensionless const& x);
}  // namespace quantities
}  // namespace principia

#include "quantities/ElementaryFunctions-body.hpp"
