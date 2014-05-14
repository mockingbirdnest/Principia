#pragma once

#include<string>

namespace principia {
namespace quantities {
// A double by any other name...
class Dimensionless {
 public:
  Dimensionless();
  // No explicit here: we want implicit conversion from double.
  Dimensionless(double const value);  // NOLINT(runtime/explicit)

  // Returns |Dimensionless(1)|, for consistency with |Quantity<D>::SIUnit()|.
  static Dimensionless SIUnit();

  double value() const;
  Dimensionless Pow(int const) const;
  // This function calls Pow(Exponent), its purpose is to provide consistency
  // with Quantity<D>.Pow<Exponent>();
  template<int Exponent>
  Dimensionless Pow() const;

 private:
  double value_;
};

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

std::string ToString(Dimensionless const& number,
                     unsigned char const precision = 16);


template<typename D>
std::ostream& operator<<(std::ostream& os, Dimensionless const& number);

}  // namespace quantities
}  // namespace principia

#include "quantities/dimensionless_body.hpp"
