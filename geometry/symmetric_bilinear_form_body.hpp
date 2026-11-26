#pragma once

#include "geometry/symmetric_bilinear_form.hpp"

#include <algorithm>
#include <limits>
#include <string>
#include <utility>

#include "numerics/elementary_functions.hpp"

namespace principia {
namespace geometry {
namespace _symmetric_bilinear_form {
namespace internal {

using namespace principia::numerics::_elementary_functions;

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
template<quantity Scalar>
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
    θ = {.cos = 1, .sin = 0};
  }
  return θ;
};


template<quantity Scalar, typename Frame, int rank>
SymmetricBilinearForm<Scalar, MultivectorTag, Frame, Rank<rank>>::
SymmetricBilinearForm(
    R3x3Matrix<Scalar> const& matrix)
    : matrix_(matrix) {
  DCHECK_EQ(matrix_(0, 1), matrix_(1, 0));
  DCHECK_EQ(matrix_(0, 2), matrix_(2, 0));
  DCHECK_EQ(matrix_(1, 2), matrix_(2, 1));
}

template<quantity Scalar, typename Frame, int rank>
SymmetricBilinearForm<Scalar, MultivectorTag, Frame, Rank<rank>>::
SymmetricBilinearForm(
    R3x3Matrix<Scalar>&& matrix)
    : matrix_(std::move(matrix)) {
  DCHECK_EQ(matrix_(0, 1), matrix_(1, 0));
  DCHECK_EQ(matrix_(0, 2), matrix_(2, 0));
  DCHECK_EQ(matrix_(1, 2), matrix_(2, 1));
}

template<quantity Scalar, typename Frame, int rank>
SymmetricBilinearForm<double, MultivectorTag, Frame, Rank<rank>> const&
SymmetricBilinearForm<Scalar, MultivectorTag, Frame, Rank<rank>>::
InnerProductForm() {
  static SymmetricBilinearForm<double, MultivectorTag, Frame, Rank<rank>> const
      form(R3x3Matrix<double>::Identity());
  return form;
}

template<quantity Scalar, typename Frame, int rank>
SymmetricBilinearForm<Scalar, MultivectorTag, Frame, Rank<rank>>&
SymmetricBilinearForm<Scalar, MultivectorTag, Frame, Rank<rank>>::
operator+=(SymmetricBilinearForm const& right) {
  return *this = *this + right;
}

template<quantity Scalar, typename Frame, int rank>
SymmetricBilinearForm<Scalar, MultivectorTag, Frame, Rank<rank>>&
SymmetricBilinearForm<Scalar, MultivectorTag, Frame, Rank<rank>>::
operator-=(SymmetricBilinearForm const& right) {
  return *this = *this - right;
}

template<quantity Scalar, typename Frame, int rank>
SymmetricBilinearForm<Scalar, MultivectorTag, Frame, Rank<rank>>&
SymmetricBilinearForm<Scalar, MultivectorTag, Frame, Rank<rank>>::
operator*=(double const right) {
  return *this = *this * right;
}

template<quantity Scalar, typename Frame, int rank>
SymmetricBilinearForm<Scalar, MultivectorTag, Frame, Rank<rank>>&
SymmetricBilinearForm<Scalar, MultivectorTag, Frame, Rank<rank>>::
operator/=(double const right) {
  return *this = *this / right;
}

template<quantity Scalar, typename Frame, int rank>
R3x3Matrix<Scalar> const&
SymmetricBilinearForm<Scalar, MultivectorTag, Frame, Rank<rank>>::
coordinates() const {
  return matrix_;
}

template<quantity Scalar, typename Frame, int rank>
template<quantity LScalar, quantity RScalar>
Product<Scalar, Product<LScalar, RScalar>>
SymmetricBilinearForm<Scalar, MultivectorTag, Frame, Rank<rank>>::
operator()(Multivector<LScalar, Frame, rank> const& left,
           Multivector<RScalar, Frame, rank> const& right) const {
  return InnerProduct(left, (*this)(_, right));
}

template<quantity Scalar, typename Frame, int rank>
template<quantity LScalar>
Multivector<Product<Scalar, LScalar>, Frame, rank>
SymmetricBilinearForm<Scalar, MultivectorTag, Frame, Rank<rank>>::operator()(
    Multivector<LScalar, Frame, rank> const& left,
    MathematicalParameterPlaceholder) const {
  return Multivector<Product<Scalar, LScalar>, Frame, rank>(
      matrix_ * left.coordinates());
}

template<quantity Scalar, typename Frame, int rank>
template<quantity RScalar>
Multivector<Product<Scalar, RScalar>, Frame, rank>
SymmetricBilinearForm<Scalar, MultivectorTag, Frame, Rank<rank>>::operator()(
    MathematicalParameterPlaceholder,
    Multivector<RScalar, Frame, rank> const& right) const {
  return Multivector<Product<Scalar, RScalar>, Frame, rank>(
      matrix_ * right.coordinates());
}

template<quantity Scalar, typename Frame, int rank>
SymmetricBilinearForm<Scalar, MultivectorTag, Frame, Rank<2>>
SymmetricBilinearForm<Scalar, MultivectorTag, Frame, Rank<rank>>::
Anticommutator() const requires(rank == 1) {
  return SymmetricBilinearForm<Scalar, MultivectorTag, Frame, Rank<2>>(
      matrix_.Trace() * R3x3Matrix<double>::Identity() - matrix_);
}

template<quantity Scalar, typename Frame, int rank>
SymmetricBilinearForm<Scalar, MultivectorTag, Frame, Rank<1>>
SymmetricBilinearForm<Scalar, MultivectorTag, Frame, Rank<rank>>::
AnticommutatorInverse() const requires(rank == 2) {
  return SymmetricBilinearForm<Scalar, MultivectorTag, Frame, Rank<1>>(
      -matrix_ + matrix_.Trace() * R3x3Matrix<double>::Identity() / 2);
}

template<quantity Scalar, typename Frame, int rank>
template<typename Eigenframe>
typename
SymmetricBilinearForm<Scalar, MultivectorTag, Frame, Rank<rank>>::
template Eigensystem<Eigenframe>
SymmetricBilinearForm<Scalar, MultivectorTag, Frame, Rank<rank>>::
Diagonalize() const {
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
      LOG(WARNING) << "Difficult diagonalization: " << matrix_
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

  return {SymmetricBilinearForm<Scalar, MultivectorTag, Eigenframe, Rank<rank>>(
              R3x3Matrix<Scalar>::DiagonalMatrix({A(0, 0), A(1, 1), A(2, 2)})),
          Rotation<Eigenframe, Frame>(v[0], v[1], v[2])};
}

template<quantity Scalar, typename Frame, int rank>
bool SymmetricBilinearForm<Scalar, MultivectorTag, Frame, Rank<rank>>::
IsDiagonalized() const {
  return matrix_.IsDiagonalMatrix() &&
         matrix_(0, 0) <= matrix_(1, 1) &&
         matrix_(1, 1) <= matrix_(2, 2);
}

template<quantity Scalar, typename Frame, int rank>
Multivector<double, Frame, rank>
SymmetricBilinearForm<Scalar, MultivectorTag, Frame, Rank<rank>>::
PrincipalAxis() const {
  using Eigenframe = geometry::_frame::Frame<struct EigenframeTag>;
  auto const eigensystem = this->template Diagonalize<Eigenframe>();
  // The last eigenvalue is the largest one.
  return eigensystem.rotation(Multivector<double, Eigenframe, rank>({0, 0, 1}));
}

template<quantity Scalar, typename Frame, int rank>
void SymmetricBilinearForm<Scalar, MultivectorTag, Frame, Rank<rank>>::
WriteToMessage(not_null<serialization::SymmetricBilinearForm*> message) const {
  Frame::WriteToMessage(message->mutable_frame());
  matrix_.WriteToMessage(message->mutable_matrix());
}

template<quantity Scalar, typename Frame, int rank>
SymmetricBilinearForm<Scalar, MultivectorTag, Frame, Rank<rank>>
SymmetricBilinearForm<Scalar, MultivectorTag, Frame, Rank<rank>>::
ReadFromMessage(serialization::SymmetricBilinearForm const& message)
  requires serializable<Frame>{
  Frame::ReadFromMessage(message.frame());
  return SymmetricBilinearForm(
      R3x3Matrix<Scalar>::ReadFromMessage(message.matrix()));
}


template<quantity Scalar>
SymmetricBilinearForm<Scalar, QuantityTag>::SymmetricBilinearForm(
    Scalar const& coordinates) : coordinates_(coordinates) {}

template<quantity Scalar>
SymmetricBilinearForm<double, QuantityTag> const&
SymmetricBilinearForm<Scalar, QuantityTag>::InnerProductForm() {
  static SymmetricBilinearForm<double, Quotient> const form(1);
  return form;
}

template<quantity Scalar>
SymmetricBilinearForm<Scalar, QuantityTag>&
SymmetricBilinearForm<Scalar, QuantityTag>::operator+=(
    SymmetricBilinearForm const& right) {
  return *this = *this + right;
}

template<quantity Scalar>
SymmetricBilinearForm<Scalar, QuantityTag>&
SymmetricBilinearForm<Scalar, QuantityTag>::operator-=(
    SymmetricBilinearForm const& right) {
  return *this = *this - right;
}

template<quantity Scalar>
SymmetricBilinearForm<Scalar, QuantityTag>&
SymmetricBilinearForm<Scalar, QuantityTag>::operator*=(
    double const right) {
  return *this = *this * right;
}

template<quantity Scalar>
SymmetricBilinearForm<Scalar, QuantityTag>&
SymmetricBilinearForm<Scalar, QuantityTag>::operator/=(
    double const right) {
  return *this = *this / right;
}

template<quantity Scalar>
Scalar const&
SymmetricBilinearForm<Scalar, QuantityTag>::coordinates() const {
  return coordinates_;
}

template<quantity Scalar>
template<quantity LScalar, quantity RScalar>
Product<Scalar, Product<LScalar, RScalar>>
SymmetricBilinearForm<Scalar, QuantityTag>::operator()(
    LScalar const& left,
    RScalar const& right) const {
  return coordinates_ * left * right;
}

template<quantity Scalar>
template<quantity LScalar>
Product<Scalar, LScalar> SymmetricBilinearForm<Scalar, QuantityTag>::operator()(
    LScalar const& left,
    MathematicalParameterPlaceholder) const {
  return coordinates_ * left;
}

template<quantity Scalar>
template<quantity RScalar>
Product<Scalar, RScalar> SymmetricBilinearForm<Scalar, QuantityTag>::operator()(
    MathematicalParameterPlaceholder,
    RScalar const& right) const {
  return coordinates_ * right;
}

template<quantity Scalar>
SymmetricBilinearForm<Scalar, QuantityTag>
SymmetricBilinearForm<Scalar, QuantityTag>::Diagonalize() const {
  return *this;
}

template<quantity Scalar>
bool SymmetricBilinearForm<Scalar, QuantityTag>::IsDiagonalized() const {
  return true;
}

template<quantity Scalar>
double SymmetricBilinearForm<Scalar, QuantityTag>::PrincipalAxis() const {
  return 1.0;
}


template<quantity Scalar, typename Tag, typename... Parameters>
SymmetricBilinearForm<Scalar, Tag, Parameters...> operator+(
    SymmetricBilinearForm<Scalar, Tag, Parameters...> const& right) {
  return SymmetricBilinearForm<Scalar, Tag, Parameters...>(
      right.coordinates());
}

template<quantity Scalar, typename Tag, typename... Parameters>
SymmetricBilinearForm<Scalar, Tag, Parameters...> operator-(
    SymmetricBilinearForm<Scalar, Tag, Parameters...> const& right) {
  return SymmetricBilinearForm<Scalar, Tag, Parameters...>(
      -right.coordinates());
}

template<quantity Scalar, typename Tag, typename... Parameters>
SymmetricBilinearForm<Scalar, Tag, Parameters...> operator+(
    SymmetricBilinearForm<Scalar, Tag, Parameters...> const& left,
    SymmetricBilinearForm<Scalar, Tag, Parameters...> const& right) {
  return SymmetricBilinearForm<Scalar, Tag, Parameters...>(
      left.coordinates() + right.coordinates());
}

template<quantity Scalar, typename Tag, typename... Parameters>
SymmetricBilinearForm<Scalar, Tag, Parameters...> operator-(
    SymmetricBilinearForm<Scalar, Tag, Parameters...> const& left,
    SymmetricBilinearForm<Scalar, Tag, Parameters...> const& right) {
  return SymmetricBilinearForm<Scalar, Tag, Parameters...>(
      left.coordinates() - right.coordinates());
}

template<convertible_to_quantity LScalar,
         quantity RScalar,
         typename Tag, typename... Parameters>
SymmetricBilinearForm<Product<LScalar, RScalar>, Tag, Parameters...>
operator*(LScalar const left,
          SymmetricBilinearForm<RScalar, Tag, Parameters...> const& right) {
  return SymmetricBilinearForm<Product<LScalar, RScalar>,
                               Tag,
                               Parameters...>(left * right.coordinates());
}

template<quantity LScalar,
         convertible_to_quantity RScalar,
         typename Tag, typename... Parameters>
SymmetricBilinearForm<Product<LScalar, RScalar>, Tag, Parameters...>
operator*(SymmetricBilinearForm<LScalar, Tag, Parameters...> const& left,
          RScalar const right) {
  return SymmetricBilinearForm<Product<LScalar, RScalar>,
                               Tag,
                               Parameters...>(left.coordinates() * right);
}

template<quantity LScalar,
         convertible_to_quantity RScalar,
         typename Tag, typename... Parameters>
SymmetricBilinearForm<Quotient<LScalar, RScalar>, Tag, Parameters...>
operator/(SymmetricBilinearForm<LScalar, Tag, Parameters...> const& left,
          RScalar const right) {
  return SymmetricBilinearForm<Quotient<LScalar, RScalar>,
                               Tag,
                               Parameters...>(left.coordinates() / right);
}

template<quantity LScalar, quantity RScalar, typename Frame, int rank>
Multivector<Quotient<LScalar, RScalar>, Frame, rank> operator/(
    Multivector<LScalar, Frame, rank> const& left,
    SymmetricBilinearForm<
        RScalar, MultivectorTag, Frame, Rank<rank>> const& right) {
  return Multivector<Quotient<LScalar, RScalar>, Frame, rank>(
      right.coordinates().Solve(left.coordinates()));
}

template<quantity LScalar, quantity RScalar, typename Frame>
SymmetricBilinearForm<Product<LScalar, RScalar>,
                      MultivectorTag, Frame, Rank<1>>
SymmetricProduct(Vector<LScalar, Frame> const& left,
                 Vector<RScalar, Frame> const& right) {
  auto const& left_coordinates = left.coordinates();
  auto const& right_coordinates = right.coordinates();
  return SymmetricBilinearForm<Product<LScalar, RScalar>,
                               MultivectorTag, Frame, Rank<1>>(
      0.5 * (KroneckerProduct(left_coordinates, right_coordinates) +
             KroneckerProduct(right_coordinates, left_coordinates)));
}

template<quantity LScalar, quantity RScalar>
SymmetricBilinearForm<Product<LScalar, RScalar>, QuantityTag>
SymmetricProduct(LScalar const& left, RScalar const& right) {
  return SymmetricBilinearForm<Product<LScalar, RScalar>, QuantityTag>(
      left * right);
}

template<quantity Scalar, typename Frame>
SymmetricBilinearForm<Square<Scalar>, MultivectorTag, Frame, Rank<1>>
SymmetricSquare(Vector<Scalar, Frame> const& vector) {
  auto const& coordinates = vector.coordinates();
  return SymmetricBilinearForm<Square<Scalar>,
                               MultivectorTag, Frame, Rank<1>>(
      KroneckerProduct(coordinates, coordinates));
}

template<quantity Scalar>
SymmetricBilinearForm<Square<Scalar>, QuantityTag>
SymmetricSquare(Scalar const& quantity) {
  return SymmetricBilinearForm<Square<Scalar>, QuantityTag>(Pow<2>(quantity));
}

template<quantity LScalar, quantity RScalar, typename Frame>
Bivector<Product<LScalar, RScalar>, Frame> Anticommutator(
    SymmetricBilinearForm<
        LScalar, MultivectorTag, Frame, Rank<1>> const& form,
    Bivector<RScalar, Frame> const& bivector) {
  auto const& form_coordinates = form.coordinates();
  auto const& bivector_coordinates = bivector.coordinates();
  return Bivector<Product<LScalar, RScalar>, Frame>(
      form_coordinates.Trace() * bivector_coordinates -
      form_coordinates * bivector_coordinates);
}

template<quantity Scalar, typename Tag, typename... Parameters>
std::string DebugString(
    SymmetricBilinearForm<Scalar, Tag, Parameters...> const& form) {
  return DebugString(form.coordinates());
}

template<quantity Scalar, typename Tag, typename... Parameters>
std::ostream& operator<<(
    std::ostream& out,
    SymmetricBilinearForm<Scalar, Tag, Parameters...> const& form) {
  out << form.coordinates();
  return out;
}

}  // namespace internal
}  // namespace _symmetric_bilinear_form
}  // namespace geometry
}  // namespace principia
