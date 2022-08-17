
#pragma once

#include "geometry/symmetric_bilinear_form.hpp"

#include <algorithm>
#include <limits>
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
inline R3x3Matrix<double> JacobiRotation(int const p,
                                         int const q,
                                         CosSin const& θ) {
  auto const& [c, s] = θ;
  R3x3Matrix<double> J = R3x3Matrix<double>::Identity();
  J(p, p) = c;
  J(q, q) = c;
  J(p, q) = s;
  J(q, p) = -s;
  return J;
};

// See [GV13] section 8.5.2, algorithm 8.5.1.
template<typename Scalar>
CosSin SymmetricSchurDecomposition2By2(R3x3Matrix<Scalar> const& A,
                                       int const p,
                                       int const q) {
  static Scalar const zero{};
  CosSin θ;
  auto& [c, s] = θ;
  if (A(p, q) != zero) {
    double const τ = (A(q, q) - A(p, p)) / (2 * A(p, q));
    double t;
    if (τ >= 0) {
      t = 1 / (τ + Sqrt(1 + τ * τ));
    } else {
      t = 1 / (τ - Sqrt(1 + τ * τ));
    }
    c = 1 / Sqrt(1 + t * t);
    s = t * c;
  } else {
    θ = {1, 0};
  }
  return θ;
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
  // when the matrix is nearly diagonal, though.
  static constexpr int max_iterations = 16;
  static constexpr double ε = std::numeric_limits<double>::epsilon() / 128;

  // [GV13], Algorithm 8.5.2.
  R3x3Matrix<Scalar> A = matrix_;
  Scalar const A_frobenius_norm = A.FrobeniusNorm();
  auto V = R3x3Matrix<double>::Identity();
  for (int k = 0; k < max_iterations; ++k) {
    Scalar max_Apq{};
    int max_p;
    int max_q;

    // Find the largest off-diagonal element and exit if it's small.
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
    if (max_Apq <= ε * A_frobenius_norm) {
      break;
    }

    auto θ = SymmetricSchurDecomposition2By2(A, max_p, max_q);
    auto const J = JacobiRotation(max_p, max_q, θ);
    A = J.Transpose() * A * J;
    V = V * J;
    if (k == max_iterations - 1) {
      LOG(ERROR) << "Difficult diagonalization: " << matrix_
                 << ", stopping with: " << A;
    }
  }

  // Now we must sort the eigenvalues by increasing value and reorder the
  // eigenvectors accordingly.  In doing so, we must track the parity of the
  // basis.
  bool odd = false;
  auto const ᵗV = V.Transpose();
  std::array<Bivector<double, Frame>, 3> v = {
      Bivector<double, Frame>(ᵗV.row_x()),
      Bivector<double, Frame>(ᵗV.row_y()),
      Bivector<double, Frame>(ᵗV.row_z())};

  auto const swap_if_needed = [&](int const i, int const j) {
    if (A(i, i) > A(j, j)) {
      odd = !odd;
      std::swap(A(i, i), A(j, j));
      std::swap(v[i], v[j]);
    }
  };

  // Unrolling a bubble sort.
  swap_if_needed(0, 1);
  swap_if_needed(1, 2);
  swap_if_needed(0, 1);

  if (odd) {
    // This choice is arbitrary but it seems to maintain pretty good
    // compatibility with the old algorithm.
    v[1] *= -1;
  }

  return {SymmetricBilinearForm<Scalar, Eigenframe, Multivector>(
              R3x3Matrix<Scalar>::DiagonalMatrix({A(0, 0), A(1, 1), A(2, 2)})),
          Rotation<Eigenframe, Frame>(v[0], v[1], v[2])};
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
