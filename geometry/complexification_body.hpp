#pragma once

#include "geometry/complexification.hpp"

namespace principia {
namespace geometry {
namespace _complexification {
namespace internal {

template<typename Vector>
template<typename V, typename>
Complexification<Vector>::Complexification(V const& real_part)
    : real_part_(real_part) {}

template<typename Vector>
Complexification<Vector>::Complexification(Vector const& real_part,
                                           Vector const& imaginary_part)
    : real_part_(real_part), imaginary_part_(imaginary_part) {}

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
typename Hilbert<Vector>::Norm²Type Complexification<Vector>::Norm²()
    const {
  return Hilbert<Vector>::Norm²(real_part_) +
         Hilbert<Vector>::Norm²(imaginary_part_);
}

template<typename Vector>
template<typename R>
Complexification<Vector>& Complexification<Vector>::operator+=(R const& right) {
  return *this = *this + right;
}

template<typename Vector>
template<typename R>
Complexification<Vector>& Complexification<Vector>::operator-=(R const& right) {
  return *this = *this - right;
}

template<typename Vector>
template<typename R>
Complexification<Vector>& Complexification<Vector>::operator*=(R const& right) {
  return *this = *this * right;
}

template<typename Vector>
template<typename R>
Complexification<Vector>& Complexification<Vector>::operator/=(R const& right) {
  return *this = *this / right;
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

template<typename L, typename Vector>
Complexification<Sum<L, Vector>> operator+(
    L const& left,
    Complexification<Vector> const& right) {
  return Complexification<Sum<L, Vector>>{left + right.real_part(),
                                          right.imaginary_part()};
}

template<typename Vector, typename R>
Complexification<Sum<Vector, R>> operator+(Complexification<Vector> const& left,
                                           R const& right) {
  return Complexification<Sum<Vector, R>>{left.real_part() + right,
                                          left.imaginary_part()};
}

template<typename Vector>
Complexification<Vector> operator+(Complexification<Vector> const& left,
                                   Complexification<Vector> const& right) {
  return Complexification<Vector>{
      left.real_part() + right.real_part(),
      left.imaginary_part() + right.imaginary_part()};
}

template<typename L, typename Vector>
Complexification<Difference<L, Vector>> operator-(
    L const& left,
    Complexification<Vector> const& right) {
  return Complexification<Difference<L, Vector>>{left - right.real_part(),
                                                 right.imaginary_part()};
}

template<typename Vector, typename R>
Complexification<Difference<Vector, R>> operator-(
    Complexification<Vector> const& left,
    R const& right) {
  return Complexification<Difference<Vector, R>>{left.real_part() - right,
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
      (lr * rr + li * ri) / denominator, (-lr * ri + li * rr) / denominator);
}

template<typename LVector, typename RVector>
Complexification<typename Hilbert<LVector, RVector>::InnerProductType>
InnerProduct(Complexification<LVector> const& left,
             Complexification<RVector> const& right) {
  return Complexification<typename Hilbert<LVector, RVector>::InnerProductType>(
      Hilbert<LVector, RVector>::InnerProduct(left.real_part(),
                                              right.real_part()) +
          Hilbert<LVector, RVector>::InnerProduct(left.imaginary_part(),
                                                  right.imaginary_part()),
      Hilbert<LVector, RVector>::InnerProduct(left.imaginary_part(),
                                              right.real_part()) -
          Hilbert<LVector, RVector>::InnerProduct(left.real_part(),
                                                  right.real_part()));
}

template<typename Vector>
std::ostream& operator<<(std::ostream& out, Complexification<Vector> const& z) {
  return out << z.real_part() << " + " << z.imaginary_part() << " i";
}

}  // namespace internal
}  // namespace _complexification
}  // namespace geometry
}  // namespace principia
