
#pragma once

#include "geometry/symmetric_bilinear_form.hpp"

#include <algorithm>
#include <string>

#include "geometry/grassmann.hpp"
#include "geometry/r3_element.hpp"
#include "geometry/rotation.hpp"
#include "geometry/sign.hpp"
#include "quantities/elementary_functions.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace geometry {
namespace internal_symmetric_bilinear_form {

using quantities::Abs;
using quantities::Angle;
using quantities::ArcCos;
using quantities::Cos;
using quantities::IsFinite;
using quantities::Sqrt;
using quantities::Square;
using quantities::si::Radian;

struct CosSin {
  double cos;
  double sin;
};

// This is J(p, q, θ) in [GV13] section 8.5.1.  This matrix is also called a
// Givens rotation.
R3x3Matrix<double> JacobiRotation(int const p,
                                  int const q,
                                  CosSin const& cos_sin) {
  R3x3Matrix<double> j = R3x3Matrix<double>::Identity();
  j(p, p) = cos_sin.cos;
  j(q, q) = cos_sin.cos;
  j(p, q) = cos_sin.sin;
  j(q, p) = -cos_sin.sin;
  return j;
};

// See [GV13] section 8.5.2.
template<typename Scalar>
CosSin SymmetricShurDecomposition2(R3x3Matrix<Scalar> const& A,
                                   int const p,
                                   int const q) {
  static Scalar const zero{};
  CosSin cos_sin;
  if (A(p, q) != zero) {
    double const τ = (A(q, q) - A(p, p)) / (2 * A(p, q));
    double t;
    if (τ >= 0) {
      t = 1 / (τ + Sqrt(1 + τ * τ));
    } else {
      t = 1 / (τ - Sqrt(1 + τ * τ));
    }
    cos_sin.cos = 1 / Sqrt(1 + t * t);
    cos_sin.sin = t * cos_sin.cos;
  } else {
    cos_sin = {1, 0};
  }
  return cos_sin;
};



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
  // As a safety measure we limit the number of iterations.  We prefer to exit
  // when the matrix is actually diagonal, though.
  static constexpr int max_iterations = 10;
  static Scalar const zero{};

  // [GV13], Algorithm 8.5.2.
  R3x3Matrix<Scalar> A = matrix_;
  auto V = R3x3Matrix<double>::Identity();
  for (int k = 0; k < max_iterations; ++k) {
    Scalar max_Apq{};
    int max_p;
    int max_q;

    // Find the largest off-diagonal element and exit if it's zero.
    for (int p = 0; p < 3; ++p) {
      for (int q = p + 1; q < 3; ++q) {
        Scalar const abs_Apq = Abs(A(p, q));
        if (abs_Apq >= max_Apq) {
          max_Apq = abs_Apq;
          max_p = p;
          max_q = q;
        }
      }
    }
    if (max_Apq == zero) {
      break;
    }

    auto cos_sin = SymmetricShurDecomposition2(A, max_p, max_q);
    auto const J = JacobiRotation(max_p, max_q, cos_sin);
    A = J.Transpose() * A * J;
    V = V * J;
    LOG(ERROR) << A << " " << V;
    if (k == max_iterations - 1) {
      LOG(WARNING) << "Difficult diagonalization: " << matrix_
                   << ", stopping with: " << A;
    }
  }

  return {SymmetricBilinearForm<Scalar, Eigenframe, Multivector>(std::move(A)),
    Rotation<Eigenframe, Frame>::Identity()};
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

  // It's possible for a column to be identically 0.  To deal with this we
  // extract the column with the largest norm.
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
