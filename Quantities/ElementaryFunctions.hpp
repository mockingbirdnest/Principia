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

template<typename D>
SquareRoot<Quantity<D>> Sqrt(Quantity<D> const& x);
}
}

#include "ElementaryFunctions-body.hpp"
