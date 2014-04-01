#pragma once

#include "Dimensionless.hpp"
#include "Quantities.hpp"

namespace Principia {
namespace Quantities {
Dimensionless Sqrt(Dimensionless const& x);
Dimensionless Log(Dimensionless const& x);
Dimensionless Log(Dimensionless const& base, Dimensionless const& x);
Dimensionless Exp(Dimensionless const& x);

Dimensionless Sin(Angle const& α);
Dimensionless Cos(Angle const& α);
Dimensionless Tan(Angle const& α);
Dimensionless Cot(Angle const& α);

Angle ArcSin(Dimensionless const& x);
Angle ArcCos(Dimensionless const& x);
Angle ArcTan(Dimensionless const& y, Dimensionless const& x = 1);
Angle ArcCot(Dimensionless const& x, Dimensionless const& y = 1);

// We consider hyperbolic functions as dealing with quotients of arc length to
// curvature radius in the hyperbolic plane, which are angles. This explains the
// use of "arc" for inverse functions.

Dimensionless Sinh(Angle const& α);
Dimensionless Cosh(Angle const& α);
Dimensionless Tanh(Angle const& α);
Dimensionless Coth(Angle const& α);

Angle ArcSinh(Dimensionless const& x);
Angle ArcCosh(Dimensionless const& x);
Angle ArcTanh(Dimensionless const& x);
Angle ArcCoth(Dimensionless const& x);

template<typename D>
SquareRoot<Quantity<D>> Sqrt(Quantity<D> const& x);
}
}

#include "ElementaryFunctions-body.hpp"
