#pragma once

#include <string>

#include "base/algebra.hpp"
#include "base/concepts.hpp"
#include "base/not_null.hpp"
#include "base/traits.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/r3x3_matrix.hpp"
#include "geometry/rotation.hpp"
#include "quantities/concepts.hpp"
#include "serialization/geometry.pb.h"

namespace principia {
namespace geometry {
namespace _symmetric_bilinear_form {
namespace internal {

using namespace principia::base::_algebra;
using namespace principia::base::_concepts;
using namespace principia::base::_not_null;
using namespace principia::base::_traits;
using namespace principia::geometry::_grassmann;
using namespace principia::geometry::_r3x3_matrix;
using namespace principia::geometry::_rotation;
using namespace principia::quantities::_concepts;

// A symmetric bilinear form with dimensionality given by `Scalar`, acting on
// vector spaces specified by the given `Generator` and `Parameters`.
template<quantity Scalar,
         template<quantity, typename...> typename Generator,
         typename... Parameters>
class SymmetricBilinearForm;

// A symmetric bilinear form with dimensionality `Scalar`, on the given kind of
// `Multivector`, expressed in the coordinates of `Frame`.
template<quantity Scalar, typename Frame, int rank>
class SymmetricBilinearForm<Scalar, MultivectorGenerator, Frame, Rank<rank>> {
 public:
  SymmetricBilinearForm() = default;
  explicit SymmetricBilinearForm(R3x3Matrix<Scalar> const& matrix);
  explicit SymmetricBilinearForm(R3x3Matrix<Scalar>&& matrix);

  friend bool operator==(SymmetricBilinearForm const& left,
                         SymmetricBilinearForm const& right) = default;
  friend bool operator!=(SymmetricBilinearForm const& left,
                         SymmetricBilinearForm const& right) = default;

  SymmetricBilinearForm& operator+=(SymmetricBilinearForm const& right);
  SymmetricBilinearForm& operator-=(SymmetricBilinearForm const& right);
  SymmetricBilinearForm& operator*=(double right);
  SymmetricBilinearForm& operator/=(double right);

  R3x3Matrix<Scalar> const& coordinates() const;

  template<quantity LScalar, quantity RScalar>
  Product<Scalar, Product<LScalar, RScalar>> operator()(
      Multivector<LScalar, Frame, rank> const& left,
      Multivector<RScalar, Frame, rank> const& right) const;

  // For a form on vectors, `Anticommutator` returns the form on bivectors
  // resulting from the commutator, i.e., up to roundoff,
  //   F.Anticommutator() * α = Anticommutator(F, α),
  //   F.Anticommutator()(α, β) = InnerProduct(α, Anticommutator(F, β)).
  // Further, note that
  //   F.Anticommutator() * Wedge(v, w) = Wedge(F * v, w) + Wedge(v, F * w),
  // which is the generalization to nonsymmetric F.
  // This operation is linear in `*this`.
  SymmetricBilinearForm<Scalar, MultivectorGenerator, Frame, Rank<2>>
  Anticommutator() const requires(rank == 1);

  // This function is the inverse of `Anticommutator()`.  It is well-defined
  // only in dimension 3, where dim ⋀²V = dim V.
  SymmetricBilinearForm<Scalar, MultivectorGenerator, Frame, Rank<1>>
  AnticommutatorInverse() const requires(rank == 2);

  // The eigensystem for a form is described by (1) the form in its eigenbasis,
  // which gives the eigenvalues; and (2) a rotation from the current basis to
  // the eigenbasis, which gives the eigenvectors.  The eigenvalues are in
  // increasing order.
  template<typename Eigenframe>
  struct Eigensystem {
    SymmetricBilinearForm<Scalar,
                          MultivectorGenerator, Eigenframe, Rank<rank>> form;
    Rotation<Eigenframe, Frame> rotation;
  };

  // Computes a form equivalent to the current one but diagonalized with
  // increasing eigenvalues.
  template<typename Eigenframe>
  Eigensystem<Eigenframe> Diagonalize() const;

  // Returns true if the form is the result of a diagonalization: its matrix is
  // diagonal and the diagonal elements are increasing.
  bool IsDiagonalized() const;

  // Returns a normalized multivector in the direction of the principal axis of
  // the form, i.e., the eigenvector with the largest eigenvalue.  If this
  // operation is executed many times, it is more efficient to call
  // `Diagonalize` and use the resulting `Eigensystem`.
  Multivector<double, Frame> PrincipalAxis() const;

  void WriteToMessage(
      not_null<serialization::SymmetricBilinearForm*> message) const;
  static SymmetricBilinearForm ReadFromMessage(
      serialization::SymmetricBilinearForm const& message)
    requires serializable<Frame>;

 private:
  // All the operations on this class must ensure that this matrix remains
  // symmetric.
  R3x3Matrix<Scalar> matrix_;

  template<quantity S,
           template<quantity, typename...> typename G,
           typename... P>
  friend class SymmetricBilinearForm;

  template<template<quantity, typename...> typename G,
           typename... P>
  friend SymmetricBilinearForm<double, G, P...> const& InnerProductForm();

  template<quantity S,
           template<quantity, typename...> typename G,
           typename... P>
  friend SymmetricBilinearForm<S, G, P...> operator+(
      SymmetricBilinearForm<S, G, P...> const& right);
  template<quantity S,
           template<quantity, typename...> typename G,
           typename... P>
  friend SymmetricBilinearForm<S, G, P...> operator-(
      SymmetricBilinearForm<S, G, P...> const& right);

  template<quantity S,
           template<quantity, typename...> typename G,
           typename... P>
  friend SymmetricBilinearForm<S, G, P...> operator+(
      SymmetricBilinearForm<S, G, P...> const& left,
      SymmetricBilinearForm<S, G, P...> const& right);
  template<quantity S,
           template<quantity, typename...> typename G,
           typename... P>
  friend SymmetricBilinearForm<S, G, P...> operator-(
      SymmetricBilinearForm<S, G, P...> const& left,
      SymmetricBilinearForm<S, G, P...> const& right);

  template<convertible_to_quantity L,
           quantity R,
           template<quantity, typename...> typename G,
           typename... P>
  friend SymmetricBilinearForm<Product<L, R>, G, P...> operator*(
      L left,
      SymmetricBilinearForm<R, G, P...> const& right);
  template<quantity L,
           convertible_to_quantity R,
           template<quantity, typename...> typename G,
           typename... P>
  friend SymmetricBilinearForm<Product<L, R>, G, P...> operator*(
      SymmetricBilinearForm<L, G, P...> const& left,
      R right);
  template<quantity L,
           convertible_to_quantity R,
           template<quantity, typename...> typename G,
           typename... P>
  friend SymmetricBilinearForm<Quotient<L, R>, G, P...> operator/(
      SymmetricBilinearForm<L, G, P...> const& left,
      R right);

  template<quantity L, quantity R, typename F, typename Rank, int r>
    requires(r == Rank::rank)
  friend _grassmann::Multivector<Product<L, R>, F, r> operator*(
      SymmetricBilinearForm<L, MultivectorGenerator, F, Rank> const& left,
      _grassmann::Multivector<R, F, r> const& right);

  template<quantity L, quantity R, typename F, typename Rank, int r>
    requires(r == Rank::rank)
  friend _grassmann::Multivector<Product<L, R>, F, r> operator*(
      _grassmann::Multivector<L, F, r> const& left,
      SymmetricBilinearForm<R, MultivectorGenerator, F, Rank> const& right);

  template<quantity L, quantity R, typename F, typename Rank, int r>
    requires(r == Rank::rank)
  friend _grassmann::Multivector<Quotient<L, R>, F, r> operator/(
      _grassmann::Multivector<L, F, r> const& left,
      SymmetricBilinearForm<R, MultivectorGenerator, F, Rank> const& right);

  template<quantity L, quantity R, typename F>
  friend SymmetricBilinearForm<Product<L, R>, MultivectorGenerator, F, Rank<1>>
  SymmetricProduct(Vector<L, F> const& left,
                   Vector<R, F> const& right);

  template<quantity L, quantity R, typename F>
  friend Bivector<Product<L, R>, F> Anticommutator(
      SymmetricBilinearForm<L, MultivectorGenerator, F, Rank<1>> const& form,
      Bivector<R, F> const& bivector);

  template<quantity S,
           template<quantity, typename...> typename G,
           typename... P>
  friend std::string DebugString(SymmetricBilinearForm<S, G, P...> const& form);

  template<quantity S,
           template<quantity, typename...> typename G,
           typename... P>
  friend std::ostream& operator<<(
      std::ostream& out,
      SymmetricBilinearForm<S, G, P...> const& form);
};

// `InnerProductForm()` is the symmetric bilinear form such that for all v and
// w, `InnerProductForm()(v, w) == InnerProduct(v, w)`.
template<template<quantity, typename...> typename Generator,
         typename... Parameters>
SymmetricBilinearForm<double, Generator, Parameters...> const&
InnerProductForm();

template<quantity Scalar,
         template<quantity, typename...> typename Generator,
         typename... Parameters>
SymmetricBilinearForm<Scalar, Generator, Parameters...> operator+(
    SymmetricBilinearForm<Scalar, Generator, Parameters...> const& right);
template<quantity Scalar,
         template<quantity, typename...> typename Generator,
         typename... Parameters>
SymmetricBilinearForm<Scalar, Generator, Parameters...> operator-(
    SymmetricBilinearForm<Scalar, Generator, Parameters...> const& right);

template<quantity Scalar,
         template<quantity, typename...> typename Generator,
         typename... Parameters>
SymmetricBilinearForm<Scalar, Generator, Parameters...> operator+(
    SymmetricBilinearForm<Scalar, Generator, Parameters...> const& left,
    SymmetricBilinearForm<Scalar, Generator, Parameters...> const& right);
template<quantity Scalar,
         template<quantity, typename...> typename Generator,
         typename... Parameters>
SymmetricBilinearForm<Scalar, Generator, Parameters...> operator-(
    SymmetricBilinearForm<Scalar, Generator, Parameters...> const& left,
    SymmetricBilinearForm<Scalar, Generator, Parameters...> const& right);

template<convertible_to_quantity LScalar,
         quantity RScalar,
         template<quantity, typename...> typename Generator,
         typename... Parameters>
SymmetricBilinearForm<Product<LScalar, RScalar>, Generator, Parameters...>
operator*(LScalar left,
          SymmetricBilinearForm<RScalar,
                                Generator, Parameters...> const& right);
template<quantity LScalar,
         convertible_to_quantity RScalar,
         template<quantity, typename...> typename Generator,
         typename... Parameters>
SymmetricBilinearForm<Product<LScalar, RScalar>, Generator, Parameters...>
operator*(SymmetricBilinearForm<LScalar, Generator, Parameters...> const& left,
          RScalar right);
template<quantity LScalar,
         convertible_to_quantity RScalar,
         template<quantity, typename...> typename Generator,
         typename... Parameters>
SymmetricBilinearForm<Quotient<LScalar, RScalar>, Generator, Parameters...>
operator/(SymmetricBilinearForm<LScalar, Generator, Parameters...> const& left,
          RScalar right);


template<quantity LScalar, quantity RScalar, typename Frame, int rank>
Multivector<Product<LScalar, RScalar>, Frame, rank> operator*(
    SymmetricBilinearForm<
        LScalar, MultivectorGenerator, Frame, Rank<rank>> const& left,
    Multivector<RScalar, Frame, rank> const& right);

template<quantity LScalar, quantity RScalar, typename Frame, int rank>
Multivector<Product<LScalar, RScalar>, Frame, rank> operator*(
    Multivector<LScalar, Frame, rank> const& left,
    SymmetricBilinearForm<
        RScalar, MultivectorGenerator, Frame, Rank<rank>> const& right);

// Solves the system `result * right == left`.  Note that by symmetry, this is
// also the solution of `right * result == left`.
template<quantity LScalar, quantity RScalar, typename Frame, int rank>
Multivector<Quotient<LScalar, RScalar>, Frame, rank> operator/(
    Multivector<LScalar, Frame, rank> const& left,
    SymmetricBilinearForm<
        RScalar, MultivectorGenerator, Frame, Rank<rank>> const& right);

// `SymmetricProduct(v, w)` is v ⊙ w ≔ (v ⊗ w + w ⊗ v) / 2.
template<quantity LScalar, quantity RScalar, typename Frame>
SymmetricBilinearForm<Product<LScalar, RScalar>,
                      MultivectorGenerator, Frame, Rank<1>>
SymmetricProduct(Vector<LScalar, Frame> const& left,
                 Vector<RScalar, Frame> const& right);

// `SymmetricSquare(v)` is `SymmetricProduct(v, v)`.
template<quantity Scalar, typename Frame>
SymmetricBilinearForm<Square<Scalar>, MultivectorGenerator, Frame, Rank<1>>
SymmetricSquare(Vector<Scalar, Frame> const& vector);

// Symmetric bilinear forms on vectors act on bivectors through this function.
// `Anticommutator(F, B)` is (tr(F)𝟙 - F)B in ℝ³ representation.  In matrix
// representation it is FB + BF = {F, B}.
template<quantity LScalar, quantity RScalar, typename Frame>
Bivector<Product<LScalar, RScalar>, Frame> Anticommutator(
    SymmetricBilinearForm<LScalar,
                          MultivectorGenerator, Frame, Rank<1>> const& form,
    Bivector<RScalar, Frame> const& bivector);

template<quantity Scalar,
         template<quantity, typename...> typename Generator,
         typename... Parameters>
std::string DebugString(
    SymmetricBilinearForm<Scalar, Generator, Parameters...> const& form);

template<quantity Scalar,
         template<quantity, typename...> typename Generator,
         typename... Parameters>
std::ostream& operator<<(
    std::ostream& out,
    SymmetricBilinearForm<Scalar, Generator, Parameters...> const& form);

}  // namespace internal

using internal::Anticommutator;
using internal::InnerProductForm;
using internal::SymmetricBilinearForm;
using internal::SymmetricProduct;
using internal::SymmetricSquare;

}  // namespace _symmetric_bilinear_form
}  // namespace geometry
}  // namespace principia

#include "geometry/symmetric_bilinear_form_body.hpp"
