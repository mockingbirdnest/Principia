
#pragma once

#include "geometry/symmetric_bilinear_form.hpp"

#include <algorithm>
#include <string>

#include "geometry/grassmann.hpp"
#include "geometry/r3_element.hpp"
#include "geometry/rotation.hpp"
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

template<typename Scalar,
         typename Frame,
         template<typename, typename> typename Multivector>
SymmetricBilinearForm<Scalar, Frame, Multivector>::SymmetricBilinearForm(
    R3x3Matrix<Scalar> const& matrix)
    : matrix_(matrix) {
  DCHECK_EQ(matrix_(0, 1), matrix_(1, 0));
  DCHECK_EQ(matrix_(0, 2), matrix_(2, 0));
  DCHECK_EQ(matrix_(1, 2), matrix_(2, 1));
}

template<typename Scalar,
         typename Frame,
         template<typename, typename> typename Multivector>
SymmetricBilinearForm<Scalar, Frame, Multivector>::SymmetricBilinearForm(
    R3x3Matrix<Scalar>&& matrix)
    : matrix_(std::move(matrix)) {
  DCHECK_EQ(matrix_(0, 1), matrix_(1, 0));
  DCHECK_EQ(matrix_(0, 2), matrix_(2, 0));
  DCHECK_EQ(matrix_(1, 2), matrix_(2, 1));
}

template<typename Scalar,
         typename Frame,
         template<typename, typename> typename Multivector>
SymmetricBilinearForm<Scalar, Frame, Multivector>&
SymmetricBilinearForm<Scalar, Frame, Multivector>::operator+=(
    SymmetricBilinearForm const& right) {
  return *this = *this + right;
}

template<typename Scalar,
         typename Frame,
         template<typename, typename> typename Multivector>
SymmetricBilinearForm<Scalar, Frame, Multivector>&
SymmetricBilinearForm<Scalar, Frame, Multivector>::operator-=(
    SymmetricBilinearForm const& right) {
  return *this = *this - right;
}

template<typename Scalar,
         typename Frame,
         template<typename, typename> typename Multivector>
SymmetricBilinearForm<Scalar, Frame, Multivector>&
SymmetricBilinearForm<Scalar, Frame, Multivector>::operator*=(
    double const right) {
  return *this = *this * right;
}

template<typename Scalar,
         typename Frame,
         template<typename, typename> typename Multivector>
SymmetricBilinearForm<Scalar, Frame, Multivector>&
SymmetricBilinearForm<Scalar, Frame, Multivector>::operator/=(
    double const right) {
  return *this = *this / right;
}

template<typename Scalar,
         typename Frame,
         template<typename, typename> typename Multivector>
R3x3Matrix<Scalar> const&
SymmetricBilinearForm<Scalar, Frame, Multivector>::coordinates()
    const {
  return matrix_;
}

template<typename Scalar,
         typename Frame,
         template<typename, typename> typename Multivector>
template<typename LScalar, typename RScalar>
Product<Scalar, Product<LScalar, RScalar>>
SymmetricBilinearForm<Scalar, Frame, Multivector>::operator()(
    Multivector<LScalar, Frame> const& left,
    Multivector<RScalar, Frame> const& right) const {
  return InnerProduct(left, *this * right);
}

template<typename Scalar,
         typename Frame,
         template<typename, typename> typename Multivector>
template<template<typename, typename> typename, typename>
SymmetricBilinearForm<Scalar, Frame, Bivector>
SymmetricBilinearForm<Scalar, Frame, Multivector>::Anticommutator() const {
  return SymmetricBilinearForm<Scalar, Frame, Bivector>(
      matrix_.Trace() * R3x3Matrix<double>::Identity() - matrix_);
}

template<typename Scalar,
         typename Frame,
         template<typename, typename> typename Multivector>
template<template<typename, typename> typename, typename>
SymmetricBilinearForm<Scalar, Frame, Vector>
SymmetricBilinearForm<Scalar, Frame, Multivector>::AnticommutatorInverse()
    const {
  return SymmetricBilinearForm<Scalar, Frame, Vector>(
      -matrix_ + matrix_.Trace() * R3x3Matrix<double>::Identity() / 2);
}

template<typename Scalar,
         typename Frame,
         template<typename, typename> typename Multivector>
template<typename Eigenframe>
typename SymmetricBilinearForm<Scalar, Frame, Multivector>::
    template Eigensystem<Eigenframe>
    SymmetricBilinearForm<Scalar, Frame, Multivector>::Diagonalize() const {
  Scalar const zero;
  R3x3Matrix<Scalar> const& A = matrix_;
  auto const I = R3x3Matrix<double>::Identity();

  // This algorithm follows
  // https://en.wikipedia.org/wiki/Eigenvalue_algorithm#3%C3%973_matrices which
  // gives closed-form formulæ for 3x3 matrices.
  Scalar const q = A.Trace() / 3;
  R3x3Matrix<Scalar> const A_minus_qI = A - q * I;
  Scalar const p = Sqrt((A_minus_qI * A_minus_qI).Trace() / 6);

  if (p == zero) {
    // A is very close to q * I.
    SymmetricBilinearForm<Scalar, Eigenframe, Multivector> form(
        R3x3Matrix<Scalar>({q, zero, zero},
                           {zero, q, zero},
                           {zero, zero, q}));
    return {form, Rotation<Eigenframe, Frame>::Identity()};
  }

  R3x3Matrix<double> const B = A_minus_qI / p;
  double const det_B = B.Determinant();
  Angle const θ = ArcCos(std::clamp(det_B, -2.0, 2.0) * 0.5);
  double const β₀ = 2 * Cos(θ / 3);
  double const β₁ = 2 * Cos((θ + 2 * π * Radian) / 3);
  double const β₂ = 2 * Cos((θ + 4 * π * Radian) / 3);
  std::array<Scalar, 3> αs = {p * β₀ + q, p * β₁ + q, p * β₂ + q};
  // We expect αs[1] <= αs[2] <= αs[0] here, but sorting ensures that we are
  // correct irrespective of numerical errors.
  std::sort(αs.begin(), αs.end());
  Scalar const& α₀ = αs[0];
  Scalar const& α₁ = αs[1];
  Scalar const& α₂ = αs[2];

  // The form in its diagonal basis.
  SymmetricBilinearForm<Scalar, Eigenframe, Multivector> form(
      R3x3Matrix<Scalar>({α₀, zero, zero},
                         {zero, α₁, zero},
                         {zero, zero, α₂}));

  // Use the Cayley-Hamilton theorem to efficiently find the eigenvectors.  The
  // m matrices contain, in columns, eigenvectors for the corresponding α.
  // However it's possible for a column to be identically 0.  To deal with this
  // the call to PickEigenvector extracts the column with the largest norm, and
  // we make sure that the eigenvectors remain orthonormal.
  R3x3Matrix<Scalar> const A_minus_α₀I = A - α₀ * I;
  R3x3Matrix<Scalar> const A_minus_α₁I = A - α₁ * I;
  R3x3Matrix<Scalar> const A_minus_α₂I = A - α₂ * I;
  auto const m₀ = A_minus_α₁I * A_minus_α₂I;
  auto const m₁ = A_minus_α₂I * A_minus_α₀I;
  auto const v₀ = Normalize(Vector<Square<Scalar>, Frame>(PickEigenvector(m₀)));
  auto const v₁ = Normalize(Vector<Square<Scalar>, Frame>(PickEigenvector(m₁))
                                .OrthogonalizationAgainst(v₀));

  Rotation<Eigenframe, Frame> const rotation{v₀, v₁, Wedge(v₀, v₁)};
  return {form, rotation};
}

template<typename Scalar,
         typename Frame,
         template<typename, typename> typename Multivector>
void SymmetricBilinearForm<Scalar, Frame, Multivector>::WriteToMessage(
    not_null<serialization::SymmetricBilinearForm*> message) const {
  Frame::WriteToMessage(message->mutable_frame());
  matrix_.WriteToMessage(message->mutable_matrix());
}

template<typename Scalar,
         typename Frame,
         template<typename, typename> typename Multivector>
template<typename, typename>
SymmetricBilinearForm<Scalar, Frame, Multivector>
SymmetricBilinearForm<Scalar, Frame, Multivector>::ReadFromMessage(
    serialization::SymmetricBilinearForm const& message) {
  Frame::ReadFromMessage(message.frame());
  return SymmetricBilinearForm(
      R3x3Matrix<Scalar>::ReadFromMessage(message.matrix()));
}

template<typename Scalar,
         typename Frame,
         template<typename, typename> typename Multivector>
template<typename S>
R3Element<S>
SymmetricBilinearForm<Scalar, Frame, Multivector>::PickEigenvector(
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
  return vs.back();
}

template<typename Frame, template<typename, typename> typename Multivector>
SymmetricBilinearForm<double, Frame, Multivector> const& InnerProductForm() {
  static auto const identity =
      SymmetricBilinearForm<double, Frame, Multivector>(
          R3x3Matrix<double>::Identity());
  return identity;
}

template<typename Scalar,
         typename Frame,
         template<typename, typename> typename Multivector>
SymmetricBilinearForm<Scalar, Frame, Multivector> operator+(
    SymmetricBilinearForm<Scalar, Frame, Multivector> const& right) {
  return SymmetricBilinearForm<Scalar, Frame, Multivector>(right.matrix_);
}

template<typename Scalar,
         typename Frame,
         template<typename, typename> typename Multivector>
SymmetricBilinearForm<Scalar, Frame, Multivector> operator-(
    SymmetricBilinearForm<Scalar, Frame, Multivector> const& right) {
  return SymmetricBilinearForm<Scalar, Frame, Multivector>(-right.matrix_);
}

template<typename Scalar,
         typename Frame,
         template<typename, typename> typename Multivector>
SymmetricBilinearForm<Scalar, Frame, Multivector> operator+(
    SymmetricBilinearForm<Scalar, Frame, Multivector> const& left,
    SymmetricBilinearForm<Scalar, Frame, Multivector> const& right) {
  return SymmetricBilinearForm<Scalar, Frame, Multivector>(left.matrix_ +
                                                           right.matrix_);
}

template<typename Scalar,
         typename Frame,
         template<typename, typename> typename Multivector>
SymmetricBilinearForm<Scalar, Frame, Multivector> operator-(
    SymmetricBilinearForm<Scalar, Frame, Multivector> const& left,
    SymmetricBilinearForm<Scalar, Frame, Multivector> const& right) {
  return SymmetricBilinearForm<Scalar, Frame, Multivector>(left.matrix_ -
                                                           right.matrix_);
}

template<typename LScalar,
         typename RScalar,
         typename Frame,
         template<typename, typename> typename Multivector>
SymmetricBilinearForm<Product<LScalar, RScalar>, Frame, Multivector> operator*(
    LScalar const left,
    SymmetricBilinearForm<RScalar, Frame, Multivector> const& right) {
  return SymmetricBilinearForm<Product<LScalar, RScalar>, Frame, Multivector>(
      left * right.matrix_);
}

template<typename LScalar,
         typename RScalar,
         typename Frame,
         template<typename, typename> typename Multivector>
SymmetricBilinearForm<Product<LScalar, RScalar>, Frame, Multivector> operator*(
    SymmetricBilinearForm<LScalar, Frame, Multivector> const& left,
    RScalar const right) {
  return SymmetricBilinearForm<Product<LScalar, RScalar>, Frame, Multivector>(
      left.matrix_ * right);
}

template<typename LScalar,
         typename RScalar,
         typename Frame,
         template<typename, typename> typename Multivector>
SymmetricBilinearForm<Quotient<LScalar, RScalar>, Frame, Multivector> operator/(
    SymmetricBilinearForm<LScalar, Frame, Multivector> const& left,
    RScalar const right) {
  return SymmetricBilinearForm<Quotient<LScalar, RScalar>, Frame, Multivector>(
      left.matrix_ / right);
}

template<typename LScalar,
         typename RScalar,
         typename Frame,
         template<typename, typename> typename M,
         int rank,
         typename>
Multivector<Product<LScalar, RScalar>, Frame, rank> operator*(
    SymmetricBilinearForm<LScalar, Frame, M> const& left,
    Multivector<RScalar, Frame, rank> const& right) {
  return Multivector<Product<LScalar, RScalar>, Frame, rank>(
      left.matrix_ * right.coordinates());
}

template<typename LScalar,
         typename RScalar,
         typename Frame,
         template<typename, typename> typename M,
         int rank,
         typename>
Multivector<Product<LScalar, RScalar>, Frame, rank> operator*(
    Multivector<LScalar, Frame, rank> const& left,
    SymmetricBilinearForm<RScalar, Frame, M> const& right) {
  return Multivector<Product<LScalar, RScalar>, Frame, rank>(
      left.coordinates() * right.matrix_);
}

template<typename LScalar,
         typename RScalar,
         typename Frame,
         template<typename, typename> typename M,
         int rank,
         typename>
Multivector<Quotient<LScalar, RScalar>, Frame, rank> operator/(
    Multivector<LScalar, Frame, rank> const& left,
    SymmetricBilinearForm<RScalar, Frame, M> const& right) {
  return Multivector<Quotient<LScalar, RScalar>, Frame, rank>(
      right.matrix_.Solve(left.coordinates()));
}

template<typename LScalar, typename RScalar, typename Frame>
SymmetricBilinearForm<Product<LScalar, RScalar>, Frame, Vector>
SymmetricProduct(Vector<LScalar, Frame> const& left,
                 Vector<RScalar, Frame> const& right) {
  return SymmetricBilinearForm<Product<LScalar, RScalar>, Frame, Vector>(
      0.5 * (KroneckerProduct(left.coordinates(), right.coordinates()) +
             KroneckerProduct(right.coordinates(), left.coordinates())));
}

template<typename LScalar, typename RScalar, typename Frame>
Bivector<Product<LScalar, RScalar>, Frame> Anticommutator(
    SymmetricBilinearForm<LScalar, Frame, Vector> const& form,
    Bivector<RScalar, Frame> const& bivector) {
  return Bivector<Product<LScalar, RScalar>, Frame>(
      form.matrix_.Trace() * bivector.coordinates() -
      form.matrix_ * bivector.coordinates());
}

template<typename Scalar,
         typename Frame,
         template<typename, typename> typename Multivector>
bool operator==(
    SymmetricBilinearForm<Scalar, Frame, Multivector> const& left,
    SymmetricBilinearForm<Scalar, Frame, Multivector> const& right) {
  return left.matrix_ == right.matrix_;
}

template<typename Scalar,
         typename Frame,
         template<typename, typename> typename Multivector>
bool operator!=(
    SymmetricBilinearForm<Scalar, Frame, Multivector> const& left,
    SymmetricBilinearForm<Scalar, Frame, Multivector> const& right) {
  return left.matrix_ != right.matrix_;
}

template<typename Scalar,
         typename Frame,
         template<typename, typename> typename Multivector>
std::string DebugString(
    SymmetricBilinearForm<Scalar, Frame, Multivector> const& form) {
  return DebugString(form.matrix_);
}

template<typename Scalar,
         typename Frame,
         template<typename, typename> typename Multivector>
std::ostream& operator<<(
    std::ostream& out,
    SymmetricBilinearForm<Scalar, Frame, Multivector> const& form) {
  out << form.matrix_;
  return out;
}

}  // namespace internal_symmetric_bilinear_form
}  // namespace geometry
}  // namespace principia
