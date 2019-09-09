
#pragma once

#include "geometry/symmetric_bilinear_form.hpp"

#include <algorithm>
#include <string>

#include "geometry/grassmann.hpp"
#include "geometry/r3_element.hpp"
#include "quantities/elementary_functions.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace geometry {
namespace internal_symmetric_bilinear_form {

using quantities::Angle;
using quantities::ArcCos;
using quantities::Cos;
using quantities::Sqrt;
using quantities::si::Radian;

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
typename SymmetricBilinearForm<Scalar, Frame>::Eigensystem<Eigenframe>
SymmetricBilinearForm<Scalar, Frame>::Diagonalize() const {
  R3x3Matrix<Scalar> const& A = matrix_;
  auto const I = R3x3Matrix<double>::Identity();

  // This algorithm follows
  // https://en.wikipedia.org/wiki/Eigenvalue_algorithm#3%C3%973_matrices which
  // gives closed-form formulæ for 3x3 matrices.
  Scalar const q = A.Trace() / 3;
  R3x3Matrix<Scalar> const A_minus_qI = A - q * I;
  Scalar const p = Sqrt((A_minus_qI * A_minus_qI).Trace() / 6);
  R3x3Matrix<double> const B = A_minus_qI / p;
  double const det_B = B.Determinant();
  Angle const θ = ArcCos(det_B * 0.5);
  double const β₀ = 2 * Cos(θ / 3);
  double const β₁ = 2 * Cos((θ + 2 * π * Radian) / 3);
  double const β₂ = 2 * Cos((θ + 4 * π * Radian) / 3);
  std::array<Scalar, 3> αs = {p * β₀ + q, p * β₁ + q, p * β₂ + q};
  std::sort(αs.begin(), αs.end());
  Scalar const& α₀ = αs[0];
  Scalar const& α₁ = αs[1];
  Scalar const& α₂ = αs[2];

  // The form in its diagonal basis.
  Scalar const zero;
  SymmetricBilinearForm<Scalar, Eigenframe> form(
      R3x3Matrix<Scalar>({  α₀, zero, zero},
                         {zero,   α₁, zero},
                         {zero, zero,   α₂}));

  // Use the Cayley-Hamilton theorem to efficiently find the eigenvectors.  The
  // m matrices contain, in columns, eigenvectors for the corresponding α.
  // However it's possible for a column to be identically 0.  To deal with this
  // the call to PickEigenvector extracts the column with the largest norm.
  R3x3Matrix<Scalar> const A_minus_α₀I = A - α₀ * I;
  R3x3Matrix<Scalar> const A_minus_α₁I = A - α₁ * I;
  R3x3Matrix<Scalar> const A_minus_α₂I = A - α₂ * I;
  auto const m₀ = A_minus_α₁I * A_minus_α₂I;
  auto const m₁ = A_minus_α₂I * A_minus_α₀I;
  auto const m₂ = A_minus_α₀I * A_minus_α₁I;
  auto const v₀ = PickEigenvector(m₀);
  auto const v₁ = PickEigenvector(m₁);

  // The vectors (v₀, v₁, v₂) forms an orthonormal basis.  Make sure that it is
  // direct.
  auto v₂ = PickEigenvector(m₂);
  if (R3x3Matrix(v₀, v₁, v₂).Determinant() < 0) {
    v₂ = -v₂;
  }
  Rotation<Frame, Eigenframe> const rotation{
      Vector<double, Frame>(v₀),
      Vector<double, Frame>(v₁),
      Bivector<double, Frame>(v₂)};  // This bivector is a bit dodgy.
  return {form, rotation};
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

template<typename Scalar, typename Frame>
template<typename S>
R3Element<double> SymmetricBilinearForm<Scalar, Frame>::PickEigenvector(
    R3x3Matrix<S> const& matrix) {
  static R3Element<double> const e₀{1, 0, 0};
  static R3Element<double> const e₁{0, 1, 0};
  static R3Element<double> const e₂{0, 0, 1};
  std::array<R3Element<S>, 3> vs = {matrix * e₀, matrix * e₁, matrix * e₂};
  std::sort(vs.begin(),
            vs.end(),
            [](R3Element<S> const& left, R3Element<S> const& right) {
              return left.Norm²() < right.Norm²();
            });
  return NormalizeOrZero(vs.back());
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
