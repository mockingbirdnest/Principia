// Dimensionless.h

#pragma once

#include<string>

namespace Principia {
namespace Quantities {
// A double by any other name...
class Dimensionless {
public:
  Dimensionless(double value);
  double Value() const;
private:
  double value_;
};

Dimensionless Exponentiate(Dimensionless const&, int const);

Dimensionless operator+(Dimensionless const&);
Dimensionless operator-(Dimensionless const&);
Dimensionless operator+(Dimensionless const&, Dimensionless const&);
Dimensionless operator-(Dimensionless const&, Dimensionless const&);
Dimensionless operator*(Dimensionless const&, Dimensionless const&);
Dimensionless operator/(Dimensionless const&, Dimensionless const&);

void operator+=(Dimensionless&, Dimensionless const&);
void operator-=(Dimensionless&, Dimensionless const&);
void operator*=(Dimensionless&, Dimensionless const&);
void operator/=(Dimensionless&, Dimensionless const&);

bool operator>(Dimensionless const&, Dimensionless const&);
bool operator<(Dimensionless const&, Dimensionless const&);
bool operator>=(Dimensionless const&, Dimensionless const&);
bool operator<=(Dimensionless const&, Dimensionless const&);
bool operator==(Dimensionless const&, Dimensionless const&);
bool operator!=(Dimensionless const&, Dimensionless const&);

Dimensionless Abs(Dimensionless const&);

std::wstring ToString(Dimensionless const&);
}
}

#include "Dimensionless.inl"
