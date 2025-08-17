#pragma once

#include <string>

#include "base/concepts.hpp"
#include "base/not_null.hpp"
#include "base/traits.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/r3x3_matrix.hpp"
#include "geometry/rotation.hpp"
#include "quantities/named_quantities.hpp"
#include "serialization/geometry.pb.h"

namespace principia {
namespace geometry {
namespace _symmetric_bilinear_form {
namespace internal {

using namespace principia::base::_concepts;
using namespace principia::base::_not_null;
using namespace principia::base::_traits;
using namespace principia::geometry::_grassmann;
using namespace principia::geometry::_r3x3_matrix;
using namespace principia::geometry::_rotation;
using namespace principia::quantities::_arithmetic;

// A symmetric bilinear form with dimensionality `Scalar`, on the given kind of
// `Multivector`, expressed in the coordinates of `Frame`.
template<typename Scalar,
         typename Frame,
         template<typename, typename> typename Multivector>
class SymmetricBilinearForm {
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

  template<typename LScalar, typename RScalar>
  Product<Scalar, Product<LScalar, RScalar>> operator()(
      Multivector<LScalar, Frame> const& left,
      Multivector<RScalar, Frame> const& right) const;

  // For a form on vectors, `Anticommutator` returns the form on bivectors
  // resulting from the commutator, i.e., up to roundoff,
  //   F.Anticommutator() * α = Anticommutator(F, α),
  //   F.Anticommutator()(α, β) = InnerProduct(α, Anticommutator(F, β)).
  // Further, note that
  //   F.Anticommutator() * Wedge(v, w) = Wedge(F * v, w) + Wedge(v, F * w),
  // which is the generalization to nonsymmetric F.
  // This operation is linear in `*this`.
  template<template<typename, typename> typename M = Multivector,
           typename = std::enable_if_t<is_same_template_v<M, Vector>>>
  SymmetricBilinearForm<Scalar, Frame, Bivector> Anticommutator() const;

  // This function is the inverse of `Anticommutator()`.  It is well-defined
  // only in dimension 3, where dim ⋀²V = dim V.
  template<template<typename, typename> typename M = Multivector,
           typename = std::enable_if_t<is_same_template_v<M, Bivector>>>
  SymmetricBilinearForm<Scalar, Frame, Vector> AnticommutatorInverse() const;

  // The eigensystem for a form is described by (1) the form in its eigenbasis,
  // which gives the eigenvalues; and (2) a rotation from the current basis to
  // the eigenbasis, which gives the eigenvectors.  The eigenvalues are in
  // increasing order.
  template<typename Eigenframe>
  struct Eigensystem {
    SymmetricBilinearForm<Scalar, Eigenframe, Multivector> form;
    Rotation<Eigenframe, Frame> rotation;
  };

  // Computes a form equivalent to the current one but diagonalized with
  // increasing eigenvalues.
  template<typename Eigenframe>
  Eigensystem<Eigenframe> Diagonalize() const;

  void WriteToMessage(
      not_null<serialization::SymmetricBilinearForm*> message) const;
  static SymmetricBilinearForm ReadFromMessage(
      serialization::SymmetricBilinearForm const& message)
    requires serializable<Frame>;

 private:
  // All the operations on this class must ensure that this matrix remains
  // symmetric.
  R3x3Matrix<Scalar> matrix_;

  template<typename S, typename F, template<typename, typename> typename M>
  friend class SymmetricBilinearForm;

  template<typename F, template<typename, typename> typename M>
  friend SymmetricBilinearForm<double, F, M> const& InnerProductForm();

  template<typename S, typename F, template<typename, typename> typename M>
  friend SymmetricBilinearForm<S, F, M> operator+(
      SymmetricBilinearForm<S, F, M> const& right);
  template<typename S, typename F, template<typename, typename> typename M>
  friend SymmetricBilinearForm<S, F, M> operator-(
      SymmetricBilinearForm<S, F, M> const& right);

  template<typename S, typename F, template<typename, typename> typename M>
  friend SymmetricBilinearForm<S, F, M> operator+(
      SymmetricBilinearForm<S, F, M> const& left,
      SymmetricBilinearForm<S, F, M> const& right);
  template<typename S, typename F, template<typename, typename> typename M>
  friend SymmetricBilinearForm<S, F, M> operator-(
      SymmetricBilinearForm<S, F, M> const& left,
      SymmetricBilinearForm<S, F, M> const& right);

  template<typename L,
           typename R,
           typename F,
           template<typename, typename> typename M>
  friend SymmetricBilinearForm<Product<L, R>, F, M> operator*(
      L left,
      SymmetricBilinearForm<R, F, M> const& right);
  template<typename L,
           typename R,
           typename F,
           template<typename, typename> typename M>
  friend SymmetricBilinearForm<Product<L, R>, F, M> operator*(
      SymmetricBilinearForm<L, F, M> const& left,
      R right);
  template<typename L,
           typename R,
           typename F,
           template<typename, typename> typename M>
  friend SymmetricBilinearForm<Quotient<L, R>, F, M> operator/(
      SymmetricBilinearForm<L, F, M> const& left,
      R right);

  template<typename L,
           typename R,
           typename F,
           template<typename, typename> typename M,
           int rank,
           typename>
  friend _grassmann::Multivector<Product<L, R>, F, rank> operator*(
      SymmetricBilinearForm<L, F, M> const& left,
      _grassmann::Multivector<R, F, rank> const& right);

  template<typename L,
           typename R,
           typename F,
           template<typename, typename> typename M,
           int rank,
           typename>
  friend _grassmann::Multivector<Product<L, R>, F, rank> operator*(
      _grassmann::Multivector<L, F, rank> const& left,
      SymmetricBilinearForm<R, F, M> const& right);

  template<typename L,
           typename R,
           typename F,
           template<typename, typename> typename M,
           int rank,
           typename>
  friend _grassmann::Multivector<Quotient<L, R>, F, rank> operator/(
      _grassmann::Multivector<L, F, rank> const& left,
      SymmetricBilinearForm<R, F, M> const& right);

  template<typename L, typename R, typename F>
  friend SymmetricBilinearForm<Product<L, R>, F, Vector>
  SymmetricProduct(Vector<L, F> const& left,
                   Vector<R, F> const& right);

  template<typename L, typename R, typename F>
  friend Bivector<Product<L, R>, F> Anticommutator(
      SymmetricBilinearForm<L, F, Vector> const& form,
      Bivector<R, F> const& bivector);

  template<typename S, typename F, template<typename, typename> typename M>
  friend std::string DebugString(SymmetricBilinearForm<S, F, M> const& form);

  template<typename S, typename F, template<typename, typename> typename M>
  friend std::ostream& operator<<(std::ostream& out,
                                  SymmetricBilinearForm<S, F, M> const& form);
};

// `InnerProductForm()` is the symmetric bilinear form such that for all v and
// w, `InnerProductForm()(v, w) == InnerProduct(v, w)`.
template<typename Frame, template<typename, typename> typename Multivector>
SymmetricBilinearForm<double, Frame, Multivector> const& InnerProductForm();

template<typename Scalar,
         typename Frame,
         template<typename, typename> typename Multivector>
SymmetricBilinearForm<Scalar, Frame, Multivector> operator+(
    SymmetricBilinearForm<Scalar, Frame, Multivector> const& right);
template<typename Scalar,
         typename Frame,
         template<typename, typename> typename Multivector>
SymmetricBilinearForm<Scalar, Frame, Multivector> operator-(
    SymmetricBilinearForm<Scalar, Frame, Multivector> const& right);

template<typename Scalar,
         typename Frame,
         template<typename, typename> typename Multivector>
SymmetricBilinearForm<Scalar, Frame, Multivector> operator+(
    SymmetricBilinearForm<Scalar, Frame, Multivector> const& left,
    SymmetricBilinearForm<Scalar, Frame, Multivector> const& right);
template<typename Scalar,
         typename Frame,
         template<typename, typename> typename Multivector>
SymmetricBilinearForm<Scalar, Frame, Multivector> operator-(
    SymmetricBilinearForm<Scalar, Frame, Multivector> const& left,
    SymmetricBilinearForm<Scalar, Frame, Multivector> const& right);

template<typename LScalar,
         typename RScalar,
         typename Frame,
         template<typename, typename> typename Multivector>
SymmetricBilinearForm<Product<LScalar, RScalar>, Frame, Multivector> operator*(
    LScalar left,
    SymmetricBilinearForm<RScalar, Frame, Multivector> const& right);
template<typename LScalar,
         typename RScalar,
         typename Frame,
         template<typename, typename> typename Multivector>
SymmetricBilinearForm<Product<LScalar, RScalar>, Frame, Multivector> operator*(
    SymmetricBilinearForm<LScalar, Frame, Multivector> const& left,
    RScalar right);
template<typename LScalar,
         typename RScalar,
         typename Frame,
         template<typename, typename> typename Multivector>
SymmetricBilinearForm<Quotient<LScalar, RScalar>, Frame, Multivector> operator/(
    SymmetricBilinearForm<LScalar, Frame, Multivector> const& left,
    RScalar right);


// NOTE(egg): An `operator*(SymmetricBilinearForm<L, F, M>, M<R, F>)` would fail
// to deduce M, for reasons that I do not quite understand (they seem to have to
// do with Vector not being the same thing as Multivector).  Instead we have
// this `enable_if` mess.

template<typename LScalar,
         typename RScalar,
         typename Frame,
         template<typename, typename> typename M,
         int rank,
         typename = std::enable_if_t<
            std::is_same_v<M<double, Frame>,
                           Multivector<double, Frame, rank>>>>
Multivector<Product<LScalar, RScalar>, Frame, rank> operator*(
    SymmetricBilinearForm<LScalar, Frame, M> const& left,
    Multivector<RScalar, Frame, rank> const& right);

template<typename LScalar,
         typename RScalar,
         typename Frame,
         template<typename, typename> typename M,
         int rank,
         typename = std::enable_if_t<
             std::is_same_v<M<double, Frame>,
                            Multivector<double, Frame, rank>>>>
Multivector<Product<LScalar, RScalar>, Frame, rank> operator*(
    Multivector<LScalar, Frame, rank> const& left,
    SymmetricBilinearForm<RScalar, Frame, M> const& right);

// Solves the system `result * right == left`.  Note that by symmetry, this is
// also the solution of `right * result == left`.
template<typename LScalar,
         typename RScalar,
         typename Frame,
         template<typename, typename> typename M,
         int rank,
         typename = std::enable_if_t<
             std::is_same_v<M<double, Frame>,
                            Multivector<double, Frame, rank>>>>
Multivector<Quotient<LScalar, RScalar>, Frame, rank> operator/(
    Multivector<LScalar, Frame, rank> const& left,
    SymmetricBilinearForm<RScalar, Frame, M> const& right);

// `SymmetricProduct(v, w)` is v ⊙ w ≔ (v ⊗ w + w ⊗ v) / 2.
template<typename LScalar, typename RScalar, typename Frame>
SymmetricBilinearForm<Product<LScalar, RScalar>, Frame, Vector>
SymmetricProduct(Vector<LScalar, Frame> const& left,
                 Vector<RScalar, Frame> const& right);

// `SymmetricSquare(v)` is `SymmetricProduct(v, v)`.
template<typename Scalar, typename Frame>
SymmetricBilinearForm<Square<Scalar>, Frame, Vector>
SymmetricSquare(Vector<Scalar, Frame> const& vector);

// Symmetric bilinear forms on vectors act on bivectors through this function.
// `Anticommutator(F, B)` is (tr(F)𝟙 - F)B in ℝ³ representation.  In matrix
// representation it is FB + BF = {F, B}.
template<typename LScalar, typename RScalar, typename Frame>
Bivector<Product<LScalar, RScalar>, Frame> Anticommutator(
    SymmetricBilinearForm<LScalar, Frame, Vector> const& form,
    Bivector<RScalar, Frame> const& bivector);

template<typename Scalar,
         typename Frame,
         template<typename, typename> typename Multivector>
std::string DebugString(
    SymmetricBilinearForm<Scalar, Frame, Multivector> const& form);

template<typename Scalar,
         typename Frame,
         template<typename, typename> typename Multivector>
std::ostream& operator<<(
    std::ostream& out,
    SymmetricBilinearForm<Scalar, Frame, Multivector> const& form);

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
