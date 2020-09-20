#pragma once

#include <type_traits>
#include <ostream>

#include "geometry/hilbert.hpp"
#include "quantities/named_quantities.hpp"

namespace principia {
namespace geometry {
namespace internal_complexification {

using quantities::Difference;
using quantities::Product;
using quantities::Quotient;
using quantities::Sum;

template<typename Vector>
class Complexification {
 public:
  explicit Complexification(Vector const& real_part);
  explicit Complexification(Vector const& real_part,
                            Vector const& imaginary_part);

  Vector const& real_part() const;
  Vector const& imaginary_part() const;

  Complexification Conjugate() const;

  typename Hilbert<Vector>::InnerProductType Norm²() const;

 private:
  Vector real_part_;
  Vector imaginary_part_;
};

template<typename Vector>
bool operator==(Complexification<Vector> const& left,
                Complexification<Vector> const& right);
template<typename Vector>
bool operator!=(Complexification<Vector> const& left,
                Complexification<Vector> const& right);

template<typename Vector>
Complexification<Vector> operator+(Complexification<Vector> const& right);
template<typename Vector>
Complexification<Vector> operator-(Complexification<Vector> const& right);

template<typename L, typename Vector>
Complexification<Sum<L, Vector>> operator+(
    L const& left,
    Complexification<Vector> const& right);
template<typename Vector, typename R>
Complexification<Sum<Vector, R>> operator+(
    Complexification<Vector> const& left,
    R const& right);
template<typename Vector>
Complexification<Vector> operator+(
    Complexification<Vector> const& left,
    Complexification<Vector> const& right);
template<typename L, typename Vector>
Complexification<Difference<L, Vector>> operator-(
    L const& left,
    Complexification<Vector> const& right);
template<typename Vector, typename R>
Complexification<Difference<Vector, R>> operator-(
    Complexification<Vector> const& left,
    R const& right);
template<typename Vector>
Complexification<Vector> operator-(
    Complexification<Vector> const& left,
    Complexification<Vector> const& right);

template<typename L, typename RVector>
Complexification<Product<L, RVector>> operator*(
    L const& left,
    Complexification<RVector> const& right);
template<typename LVector, typename R>
Complexification<Product<LVector, R>> operator*(
    Complexification<LVector> const& left,
    R const& right);
// TODO(egg): Numerical analysis.
template<typename LVector, typename RVector>
Complexification<Product<LVector, RVector>> operator*(
    Complexification<LVector> const& left,
    Complexification<RVector> const& right);

// TODO(egg): Numerical analysis.
template<typename L, typename RVector>
Complexification<Quotient<L, RVector>> operator/(
    L const& left,
    Complexification<RVector> const& right);
template<typename LVector, typename R>
Complexification<Quotient<LVector, R>> operator/(
    Complexification<LVector> const& left,
    R const& right);
// TODO(egg): Numerical analysis.
template<typename LVector, typename RVector>
Complexification<Quotient<LVector, RVector>> operator/(
    Complexification<LVector> const& left,
    Complexification<RVector> const& right);

template<typename Vector>
std::ostream& operator<<(std::ostream& out, Complexification<Vector> const& z);

}  // namespace internal_complexification

using internal_complexification::Complexification;

}  // namespace geometry
}  // namespace principia

#include "geometry/complexification_body.hpp"
