
#pragma once

#include "geometry/symmetric_bilinear_form.hpp"

#include <string>

namespace principia {
namespace geometry {
namespace internal_symmetric_bilinear_form {

template<typename Scalar, typename Frame>
SymmetricBilinearForm<Scalar, Frame>& SymmetricBilinearForm<Scalar, Frame>::
operator+=(SymmetricBilinearForm const& right) {
  return *this = *this + right;
}

template<typename Scalar, typename Frame>
SymmetricBilinearForm<Scalar, Frame>& SymmetricBilinearForm<Scalar, Frame>::
operator-=(SymmetricBilinearForm const& right) {
  return *this = *this - right;
}

template<typename Scalar, typename Frame>
SymmetricBilinearForm<Scalar, Frame>& SymmetricBilinearForm<Scalar, Frame>::
operator*=(double const right) {
  return *this = *this * right;
}

template<typename Scalar, typename Frame>
SymmetricBilinearForm<Scalar, Frame>& SymmetricBilinearForm<Scalar, Frame>::
operator/=(double const right) {
  return *this = *this / right;
}

template<typename Scalar, typename Frame>
template<typename LScalar, typename RScalar>
Product<Scalar, Product<LScalar, RScalar>>
SymmetricBilinearForm<Scalar, Frame>::operator()(
    Vector<LScalar, Frame> const& left,
    Vector<RScalar, Frame> const& right) const {
  return InnerProduct(left, *this * right);
}

template<typename Scalar, typename Frame>
template<typename LScalar, typename RScalar>
Product<Scalar, Product<LScalar, RScalar>>
SymmetricBilinearForm<Scalar, Frame>::operator()(
    Bivector<LScalar, Frame> const& left,
    Bivector<RScalar, Frame> const& right) const {
  return InnerProduct(left, *this * right);
}

template<typename Scalar, typename Frame>
SymmetricBilinearForm<Scalar, Frame>::SymmetricBilinearForm(
    R3x3Matrix<Scalar> const& matrix) : matrix_(matrix) {
  DCHECK_EQ(matrix_(0, 1), matrix_(1, 0));
  DCHECK_EQ(matrix_(0, 2), matrix_(2, 0));
  DCHECK_EQ(matrix_(1, 2), matrix_(2, 1));
}

template<typename Scalar, typename Frame>
SymmetricBilinearForm<Scalar, Frame>::SymmetricBilinearForm(
    R3x3Matrix<Scalar>&& matrix) : matrix_(std::move(matrix)) {
  DCHECK_EQ(matrix_(0, 1), matrix_(1, 0));
  DCHECK_EQ(matrix_(0, 2), matrix_(2, 0));
  DCHECK_EQ(matrix_(1, 2), matrix_(2, 1));
}

template<typename Frame>
SymmetricBilinearForm<double, Frame> const& InnerProductForm() {
  static auto const identity =
      SymmetricBilinearForm<double, Frame>(R3x3Matrix<double>::Identity());
  return identity;
}

template<typename Scalar, typename Frame>
SymmetricBilinearForm<Scalar, Frame> operator+(
    SymmetricBilinearForm<Scalar, Frame> const& right) {
  return SymmetricBilinearForm<Scalar, Frame>(right.matrix_);
}

template<typename Scalar, typename Frame>
SymmetricBilinearForm<Scalar, Frame> operator-(
    SymmetricBilinearForm<Scalar, Frame> const& right) {
  return SymmetricBilinearForm<Scalar, Frame>(-right.matrix_);
}

template<typename Scalar, typename Frame>
SymmetricBilinearForm<Scalar, Frame> operator+(
    SymmetricBilinearForm<Scalar, Frame> const& left,
    SymmetricBilinearForm<Scalar, Frame> const& right) {
  return SymmetricBilinearForm<Scalar, Frame>(left.matrix_ + right.matrix_);
}

template<typename Scalar, typename Frame>
SymmetricBilinearForm<Scalar, Frame> operator-(
    SymmetricBilinearForm<Scalar, Frame> const& left,
    SymmetricBilinearForm<Scalar, Frame> const& right) {
  return SymmetricBilinearForm<Scalar, Frame>(left.matrix_ - right.matrix_);
}

template<typename Scalar, typename Frame>
SymmetricBilinearForm<Scalar, Frame> operator*(
    double const left,
    SymmetricBilinearForm<Scalar, Frame> const& right) {
  return SymmetricBilinearForm<Scalar, Frame>(left * right.matrix_);
}

template<typename Scalar, typename Frame>
SymmetricBilinearForm<Scalar, Frame> operator*(
    SymmetricBilinearForm<Scalar, Frame> const& left,
    double const right) {
  return SymmetricBilinearForm<Scalar, Frame>(left.matrix_ * right);
}

template<typename Scalar, typename Frame>
SymmetricBilinearForm<Scalar, Frame> operator/(
    SymmetricBilinearForm<Scalar, Frame> const& left,
    double const right) {
  return SymmetricBilinearForm<Scalar, Frame>(left.matrix_ / right);
}

template<typename LScalar, typename RScalar, typename Frame>
Vector<Product<LScalar, RScalar>, Frame> operator*(
    SymmetricBilinearForm<LScalar, Frame> const& left,
    Vector<RScalar, Frame> const& right) {
  return Vector<Product<LScalar, RScalar>, Frame>(left.matrix_ *
                                                  right.coordinates());
}

template<typename LScalar, typename RScalar, typename Frame>
Bivector<Product<LScalar, RScalar>, Frame> operator*(
    SymmetricBilinearForm<LScalar, Frame> const& left,
    Bivector<RScalar, Frame> const& right) {
  return Bivector<Product<LScalar, RScalar>, Frame>(left.matrix_ *
                                                    right.coordinates());
}

template<typename LScalar, typename RScalar, typename Frame>
Vector<Product<LScalar, RScalar>, Frame> operator*(
    Vector<LScalar, Frame> const& left,
    SymmetricBilinearForm<RScalar, Frame> const& right) {
  return Vector<Product<LScalar, RScalar>, Frame>(left.coordinates() *
                                                  right.matrix_);
}

template<typename LScalar, typename RScalar, typename Frame>
Bivector<Product<LScalar, RScalar>, Frame> operator*(
    Bivector<LScalar, Frame> const& left,
    SymmetricBilinearForm<RScalar, Frame> const& right) {
  return Bivector<Product<LScalar, RScalar>, Frame>(left.coordinates() *
                                                    right.matrix_);
}

template<typename LScalar, typename RScalar, typename Frame>
SymmetricBilinearForm<Product<LScalar, RScalar>, Frame> SymmetricProduct(
    Vector<LScalar, Frame> const& left,
    Vector<RScalar, Frame> const& right) {
  return SymmetricBilinearForm<Product<LScalar, RScalar>, Frame>(
      0.5 * (KroneckerProduct(left.coordinates(), right.coordinates()) +
             KroneckerProduct(right.coordinates(), left.coordinates())));
}

template<typename LScalar, typename RScalar, typename Frame>
SymmetricBilinearForm<Product<LScalar, RScalar>, Frame> SymmetricProduct(
    Bivector<LScalar, Frame> const& left,
    Bivector<RScalar, Frame> const& right) {
  return SymmetricBilinearForm<Product<LScalar, RScalar>, Frame>(
      0.5 * (KroneckerProduct(left.coordinates(), right.coordinates()) +
             KroneckerProduct(right.coordinates(), left.coordinates())));
}

template<typename Scalar, typename Frame>
bool operator==(SymmetricBilinearForm<Scalar, Frame> const& left,
                SymmetricBilinearForm<Scalar, Frame> const& right) {
  return left.matrix_ == right.matrix_;
}

template<typename Scalar, typename Frame>
bool operator!=(SymmetricBilinearForm<Scalar, Frame> const& left,
                SymmetricBilinearForm<Scalar, Frame> const& right) {
  return left.matrix_ != right.matrix_;
}

template<typename Scalar, typename Frame>
std::string DebugString(SymmetricBilinearForm<Scalar, Frame> const& form) {
  return DebugString(form.matrix_);
}

template<typename Scalar, typename Frame>
std::ostream& operator<<(std::ostream& out,
                         SymmetricBilinearForm<Scalar, Frame> const& form) {
  out << form.matrix_;
  return out;
}

}  // namespace internal_symmetric_bilinear_form
}  // namespace geometry
}  // namespace principia
