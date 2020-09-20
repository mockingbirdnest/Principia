#pragma once

#include "geometry/complexification.hpp"

namespace principia {
namespace geometry {
namespace internal_complexification {

template<typename Vector>
Complexification<Vector>::Complexification(Vector const& real_part)
    : real_part_(real_part) {}

template<typename Vector>
Complexification<Vector>::Complexification(Vector const& real_part,
                                           Vector const& imaginary_part)
    : real_part_(real_part), imaginary_part_(imaginary_part) {}

template<typename Vector>
template<typename V, typename>
Complexification<Vector>::Complexification(std::complex<V> const& z)
    : real_part_(z.real()), imaginary_part_(z.imag()) {}

template<typename Vector>
Vector const& Complexification<Vector>::real_part() const {
  return real_part_;
}

template<typename Vector>
Vector const& Complexification<Vector>::imaginary_part() const {
  return imaginary_part_;
}

template<typename Vector>
Complexification<Vector> Complexification<Vector>::Conjugate() const {
  return Complexification{real_part(), -imaginary_part()};
}

template<typename Vector>
typename Hilbert<Vector>::InnerProductType Complexification<Vector>::Norm²()
    const {
  // TODO(egg): Hilbert::Norm².
  return Hilbert<Vector>::InnerProduct(real_part_, real_part_) +
         Hilbert<Vector>::InnerProduct(imaginary_part_, imaginary_part_);
}

template<typename Vector>
bool operator==(Complexification<Vector> const& left,
                Complexification<Vector> const& right) {
  return left.real_part() == right.real_part() &&
         left.imaginary_part() == right.imaginary_part();
}

template<typename Vector>
bool operator!=(Complexification<Vector> const& left,
                Complexification<Vector> const& right) {
  return left.real_part() != right.real_part() ||
         left.imaginary_part() != right.imaginary_part();
}

template<typename Vector>
Complexification<Vector> operator+(Complexification<Vector> const& right) {
  return right;
}

template<typename Vector>
Complexification<Vector> operator-(Complexification<Vector> const& right) {
  return Complexification<Vector>{-right.real_part(), -right.imaginary_part()};
}

template<typename Vector>
Complexification<Vector> operator+(Vector const& left,
                                   Complexification<Vector> const& right) {
  return Complexification<Vector>{left + right.real_part(),
                                  right.imaginary_part()};
}

template<typename Vector>
Complexification<Vector> operator+(Complexification<Vector> const& left,
                                   Vector const& right) {
  return Complexification<Vector>{left.real_part() + right,
                                  left.imaginary_part()};
}

template<typename Vector>
Complexification<Vector> operator+(Complexification<Vector> const& left,
                                   Complexification<Vector> const& right) {
  return Complexification<Vector>{
      left.real_part() + right.real_part(),
      left.imaginary_part() + right.imaginary_part()};
}

template<typename Vector>
Complexification<Vector> operator-(Vector const& left,
                                   Complexification<Vector> const& right) {
  return Complexification<Vector>{left - right.real_part(),
                                  -right.imaginary_part()};
}

template<typename Vector>
Complexification<Vector> operator-(Complexification<Vector> const& left,
                                   Vector const& right) {
  return Complexification<Vector>{left.real_part() - right,
                                  left.imaginary_part()};
}

template<typename Vector>
Complexification<Vector> operator-(Complexification<Vector> const& left,
                                   Complexification<Vector> const& right) {
  return Complexification<Vector>{
      left.real_part() - right.real_part(),
      left.imaginary_part() - right.imaginary_part()};
}

template<typename L, typename RVector>
Complexification<Product<L, RVector>> operator*(
    L const& left,
    Complexification<RVector> const& right) {
  return Complexification<Product<L, RVector>>(left * right.real_part(),
                                               left * right.imaginary_part());
}

template<typename LVector, typename R>
Complexification<Product<LVector, R>> operator*(
    Complexification<LVector> const& left,
    R const& right) {
  return Complexification<Product<LVector, R>>(left.real_part() * right,
                                               left.imaginary_part() * right);
}

template<typename LVector, typename RVector>
Complexification<Product<LVector, RVector>> operator*(
    Complexification<LVector> const& left,
    Complexification<RVector> const& right) {
  auto const& lr = left.real_part();
  auto const& li = left.imaginary_part();
  auto const& rr = right.real_part();
  auto const& ri = right.imaginary_part();
  return Complexification<Product<LVector, RVector>>(lr * rr - li * ri,
                                                     lr * ri + li * rr);
}

template<typename L, typename RVector>
Complexification<Quotient<L, RVector>> operator/(
    L const& left,
    Complexification<RVector> const& right) {
  auto const& rr = right.real_part();
  auto const& ri = right.imaginary_part();
  auto const& denominator = rr * rr + ri * ri;
  return Complexification<Quotient<L, RVector>>(left * rr / denominator,
                                                -left * ri / denominator);
}

template<typename LVector, typename R>
Complexification<Quotient<LVector, R>> operator/(
    Complexification<LVector> const& left,
    R const& right) {
  return Complexification<Quotient<LVector, R>>{left.real_part() / right,
                                                left.imaginary_part() / right};
}

template<typename LVector, typename RVector>
Complexification<Quotient<LVector, RVector>> operator/(
    Complexification<LVector> const& left,
    Complexification<RVector> const& right) {
  auto const& lr = left.real_part();
  auto const& li = left.imaginary_part();
  auto const& rr = right.real_part();
  auto const& ri = right.imaginary_part();
  auto const& denominator = rr * rr + ri * ri;
  return Complexification<Quotient<LVector, RVector>>(
      (lr * rr + li * ri) / denominator,
      (-lr * ri + li * rr) / denominator);
}

template<typename Vector>
std::ostream& operator<<(std::ostream& out, Complexification<Vector> const& z) {
  return out << z.real_part() << " + i " << z.imaginary_part();
}

}  // namespace internal_complexification
}  // namespace geometry
}  // namespace principia
