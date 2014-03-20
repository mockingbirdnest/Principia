// Constants.h

#pragma once

#include "PhysicalQuantities.h"

namespace PhysicalQuantities {

const Entropy BoltzmannConstant = 1.3806488e-23 * (Joule / Kelvin);

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
  double numberOfKelvins = Value((9.8 * Kelvin - Celsius(14)) / (1.0 * Kelvin));
  DimensionlessScalar N = Dimensionless(1e23);
  Volume V = 5 * (Metre * Metre * Metre) + 2 * Litre;
  Temperature T = 3 * Kelvin;
  Pressure P = N * BoltzmannConstant * T / V;
}

}