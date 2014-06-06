#pragma once

#include<cfloat>
#include<string>

namespace principia {
namespace quantities {
// A double by any other name...
class Dimensionless {
 public:
  typedef Dimensionless Inverse;

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

  template<typename D>
  friend class Quantity;
  template<typename D>
  friend Quantity<D> operator*(Quantity<D> const&, Dimensionless const&);
  template<typename D>
  friend Quantity<D> operator*(Dimensionless const&, Quantity<D> const&);
  template<typename D>
  friend Quantity<D> operator/(Quantity<D> const&, Dimensionless const&);
  template<typename D>
  friend typename Quantity<D>::Inverse operator/(Dimensionless const&,
                                                 Quantity<D> const&);
  template<typename D>
  friend void operator*=(Quantity<D>&, Dimensionless const&);
  template<typename D>
  friend void operator/=(Quantity<D>&, Dimensionless const&);
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
Dimensionless Max(Dimensionless const&, Dimensionless const&);

std::string ToString(Dimensionless const& number,
                     unsigned char const precision = DBL_DIG + 1);

template<typename D>
std::ostream& operator<<(std::ostream& out, Dimensionless const& number);

}  // namespace quantities
}  // namespace principia

#include "quantities/dimensionless_body.hpp"
