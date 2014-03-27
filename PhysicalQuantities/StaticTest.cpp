// StaticTest.h

#include "stdafx.h"

#include "StaticTest.h"

#include "PhysicalQuantities.h"
#include "NamedQuantities.h"
#include "SIUnits.h"
#include "Constants.h"
namespace PhysicalQuantitiesTest {
using namespace PhysicalQuantities;
public ref class staticTest {
  public:
  double test() {
    Dimensionless z = 42;
    Dimensionless y = 75;
    Dimensionless aoeu = z * y;
    Mass m = 5.0 * Kilogram;
    Speed v = 1.2 * (Metre / Second);
    v += 43 * (Metre / Second);
    v *= 5.2;
    v /= 0.7;
    Dimensionless x = 3.0;
    Dimensionless xx = 5 * x;
    Momentum p = m * v;
    Mass M = 5.0 * m;
    Energy E = Dimensionless(.5) * m * v * v;
    Force F = 1000.0 * Newton;
    double numberOfKelvins
      = ((9.8 * Kelvin - 12 * Joule / BoltzmannConstant) / (1.0 * Kelvin)).Value();
    Dimensionless N = 1e23;
    Volume V = 5 * (Metre * Metre * Metre) + 2 * Litre;
    Temperature T = 3 * Kelvin;
    Pressure P = N * BoltzmannConstant * T / V;
    MomentOfInertia I = 1 * Metre * Metre * 4 * Kilogram / (Radian * Radian);
    AngularMomentum L = I * (4 * (Radian / Second));
    return (m / (1 * Kilogram) + E / (1 * Joule) + P / (1 * Pascal)).Value()
      + numberOfKelvins;
  }
};
}