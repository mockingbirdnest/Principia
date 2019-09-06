
#pragma once

#include "geometry/symmetric_bilinear_form.hpp"

#include <string>

#include "quantities/elementary_functions.hpp"
#include "quantities/quantities.hpp"

namespace principia {
namespace geometry {
namespace internal_symmetric_bilinear_form {

using quantities::Angle;
using quantities::ArcCos;
using quantities::Cos;
using quantities::Sqrt;

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
template<typename Eigenframe>
SymmetricBilinearForm<Scalar, Frame>::Eigensystem<Eigenframe>
SymmetricBilinearForm<Scalar, Frame>::Diagonalize() const {
  // https://en.wikipedia.org/wiki/Eigenvalue_algorithm#3%C3%973_matrices
  R3x3Matrix<Scalar> const& A = matrix_;
  Scalar const q = A.Trace() / 3;
  R3x3Matrix<Scalar> const A_minus_qI = A - q * R3x3Matrix<double>::Identity();
  Scalar const p = Sqrt((A_minus_qI * A_minus_qI).Trace() / 6);
  R3x3Matrix<double> const B = A_minus_qI / p;
  Scalar const det_B = B.Determinant();
  Angle const θ = ArcCos(det_B * 0.5);
  double const β₀ = 2 * Cos(θ / 3);
  double const β₁ = 2 * Cos((θ + 2 * π) / 3);
  double const β₂ = 2 * Cos((θ + 4 * π) / 3);
}

template<typename Scalar, typename Frame>
void SymmetricBilinearForm<Scalar, Frame>::WriteToMessage(
    not_null<serialization::SymmetricBilinearForm*> message) const {
  Frame::WriteToMessage(message->mutable_frame());
  matrix_.WriteToMessage(message->mutable_matrix());
}

template<typename Scalar, typename Frame>
SymmetricBilinearForm<Scalar, Frame>
SymmetricBilinearForm<Scalar, Frame>::ReadFromMessage(
    serialization::SymmetricBilinearForm const& message) {
  Frame::ReadFromMessage(message.frame());
  return SymmetricBilinearForm(
      R3x3Matrix<Scalar>::ReadFromMessage(message.matrix()));
}

template<typename Scalar, typename Frame>
SymmetricBilinearForm<Scalar, Frame>::SymmetricBilinearForm(
    R3x3Matrix<Scalar> const& matrix)
    : matrix_(matrix) {
  DCHECK_EQ(matrix_(0, 1), matrix_(1, 0));
  DCHECK_EQ(matrix_(0, 2), matrix_(2, 0));
  DCHECK_EQ(matrix_(1, 2), matrix_(2, 1));
}

template<typename Scalar, typename Frame>
SymmetricBilinearForm<Scalar, Frame>::SymmetricBilinearForm(
    R3x3Matrix<Scalar>&& matrix)
    : matrix_(std::move(matrix)) {
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
Vector<Product<LScalar, RScalar>, Frame> operator*(
    Vector<LScalar, Frame> const& left,
    SymmetricBilinearForm<RScalar, Frame> const& right) {
  return Vector<Product<LScalar, RScalar>, Frame>(left.coordinates() *
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
Bivector<Product<LScalar, RScalar>, Frame> Anticommutator(
    SymmetricBilinearForm<LScalar, Frame> const& form,
    Bivector<RScalar, Frame> const& bivector) {
  return Bivector<Product<LScalar, RScalar>, Frame>(
      form.matrix_.Trace() * bivector.coordinates() -
      form.matrix_ * bivector.coordinates());
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
