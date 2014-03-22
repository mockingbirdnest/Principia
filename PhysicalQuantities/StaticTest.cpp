// StaticTest.h

#include "stdafx.h"

#include "StaticTest.h"

#include "PhysicalQuantities.h"
#include "NamedQuantities.h"
#include "Units.h"
#include "Constants.h"

using namespace PhysicalQuantities;

void test() {
  Mass m = 5.0 * Kilogram;
  Speed v = 1.2 * (Metre / Second);
  v += 43 * (Metre / Second);
  v *= Dimensionless(5.2);
  v /= Dimensionless(0.7);
  DimensionlessScalar x = Dimensionless(3.0);
  Momentum p = m * v;
  Energy E = Dimensionless(.5) * m * v * v;
  Force F = 1000.0 * Newton;
  double numberOfKelvins
    = Value((9.8 * Kelvin - 12 * Joule / BoltzmannConstant) / (1.0 * Kelvin));
  DimensionlessScalar N = Dimensionless(1e23);
  Volume V = 5 * (Metre * Metre * Metre) + 2 * Litre;
  Temperature T = 3 * Kelvin;
  Pressure P = N * BoltzmannConstant * T / V;
  MomentOfInertia I = (1 * (Metre * Metre)) * (4 * Kilogram);
  AngularMomentum L = I * (4 * (Radian / Second));
}