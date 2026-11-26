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
#include "quantities/quantities.hpp"
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
using namespace principia::quantities::_quantities;

// A symmetric bilinear form with dimensionality given by `Scalar`, acting on
// vector spaces parameterized by `Parameters`.  The `Tag` makes it possible to
// distinguish specializations that have the same `Parameters`.
template<quantity Scalar, typename Tag, typename... Parameters>
class SymmetricBilinearForm;

// A symmetric bilinear form with dimensionality `Scalar`, on the given kind of
// `Multivector`, expressed in the coordinates of `Frame`.
template<quantity Scalar, typename Frame, int rank>
class SymmetricBilinearForm<Scalar, MultivectorTag, Frame, Rank<rank>> {
 public:
  SymmetricBilinearForm() = default;
  explicit SymmetricBilinearForm(R3x3Matrix<Scalar> const& matrix);
  explicit SymmetricBilinearForm(R3x3Matrix<Scalar>&& matrix);

  // `InnerProductForm()` is the symmetric bilinear form such that for all v and
  // w, `InnerProductForm()(v, w) == InnerProduct(v, w)`.
  static SymmetricBilinearForm<double, MultivectorTag, Frame, Rank<rank>> const&
  InnerProductForm();

  friend bool operator==(SymmetricBilinearForm const& left,
                         SymmetricBilinearForm const& right) = default;
  friend bool operator!=(SymmetricBilinearForm const& left,
                         SymmetricBilinearForm const& right) = default;

  SymmetricBilinearForm& operator+=(SymmetricBilinearForm const& right);
  SymmetricBilinearForm& operator-=(SymmetricBilinearForm const& right);
  SymmetricBilinearForm& operator*=(double right);
  SymmetricBilinearForm& operator/=(double right);

  R3x3Matrix<Scalar> const& coordinates() const;

  // Application of the form to a pair of multivectors.
  template<quantity LScalar, quantity RScalar>
  Product<Scalar, Product<LScalar, RScalar>> operator()(
      Multivector<LScalar, Frame, rank> const& left,
      Multivector<RScalar, Frame, rank> const& right) const;

  // Application of the form to a single multivector, yielding a linear form
  // which we identify with a multivector of the same rank.
  template<quantity LScalar>
  Multivector<Product<Scalar, LScalar>, Frame, rank> operator()(
      Multivector<LScalar, Frame, rank> const& left,
      MathematicalParameterPlaceholder) const;
  template<quantity RScalar>
  Multivector<Product<Scalar, RScalar>, Frame, rank> operator()(
      MathematicalParameterPlaceholder,
      Multivector<RScalar, Frame, rank> const& right) const;

  // For a form on vectors, `Anticommutator` returns the form on bivectors
  // resulting from the commutator, i.e., up to roundoff,
  //   F.Anticommutator() * α = Anticommutator(F, α),
  //   F.Anticommutator()(α, β) = InnerProduct(α, Anticommutator(F, β)).
  // Further, note that
  //   F.Anticommutator() * Wedge(v, w) = Wedge(F * v, w) + Wedge(v, F * w),
  // which is the generalization to nonsymmetric F.
  // This operation is linear in `*this`.
  SymmetricBilinearForm<Scalar, MultivectorTag, Frame, Rank<2>>
  Anticommutator() const requires(rank == 1);

  // This function is the inverse of `Anticommutator()`.  It is well-defined
  // only in dimension 3, where dim ⋀²V = dim V.
  SymmetricBilinearForm<Scalar, MultivectorTag, Frame, Rank<1>>
  AnticommutatorInverse() const requires(rank == 2);

  // The eigensystem for a form is described by (1) the form in its eigenbasis,
  // which gives the eigenvalues; and (2) a rotation from the current basis to
  // the eigenbasis, which gives the eigenvectors.  The eigenvalues are in
  // increasing order.
  template<typename Eigenframe>
  struct Eigensystem {
    SymmetricBilinearForm<Scalar, MultivectorTag, Eigenframe, Rank<rank>> form;
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
  Multivector<double, Frame, rank> PrincipalAxis() const;

  void WriteToMessage(
      not_null<serialization::SymmetricBilinearForm*> message) const;
  static SymmetricBilinearForm ReadFromMessage(
      serialization::SymmetricBilinearForm const& message)
    requires serializable<Frame>;

 private:
  // All the operations on this class must ensure that this matrix remains
  // symmetric.
  R3x3Matrix<Scalar> matrix_;
};


// A symmetric bilinear form with dimensionality `Scalar`, on quantities.
template<quantity Scalar>
class SymmetricBilinearForm<Scalar, QuantityTag> {
 public:
  SymmetricBilinearForm() = default;
  explicit SymmetricBilinearForm(Scalar const& coordinates);

  // Returns a form with `coordinates` set to the unit of `Scalar`.
  static SymmetricBilinearForm<double, QuantityTag> const&
  InnerProductForm();

  friend bool operator==(SymmetricBilinearForm const& left,
                         SymmetricBilinearForm const& right) = default;
  friend bool operator!=(SymmetricBilinearForm const& left,
                         SymmetricBilinearForm const& right) = default;

  SymmetricBilinearForm& operator+=(SymmetricBilinearForm const& right);
  SymmetricBilinearForm& operator-=(SymmetricBilinearForm const& right);
  SymmetricBilinearForm& operator*=(double right);
  SymmetricBilinearForm& operator/=(double right);

  Scalar const& coordinates() const;

  // Application of the form to a pair of quantities (considered as members of a
  // vector space).
  template<quantity LScalar, quantity RScalar>
  Product<Scalar, Product<LScalar, RScalar>> operator()(
      LScalar const& left,
      RScalar const& right) const;

  // Application of the form to a single multivector, yielding a linear form
  // which we identify with an element of the quantities vector space.
  template<quantity LScalar>
  Product<Scalar, LScalar> operator()(LScalar const& left,
                                      MathematicalParameterPlaceholder) const;
  template<quantity RScalar>
  Product<Scalar, RScalar> operator()(MathematicalParameterPlaceholder,
                                      RScalar const& right) const;

  // Returns a copy of this object.
  SymmetricBilinearForm<Scalar, QuantityTag> Diagonalize() const;

  // Always returns true.
  bool IsDiagonalized() const;

  // Always returns 1.
  double PrincipalAxis() const;

  // No serialization support for now.

 private:
  Scalar coordinates_;
};


// Additive group.

template<quantity Scalar, typename Tag, typename... Parameters>
SymmetricBilinearForm<Scalar, Tag, Parameters...> operator+(
    SymmetricBilinearForm<Scalar, Tag, Parameters...> const& right);
template<quantity Scalar, typename Tag, typename... Parameters>
SymmetricBilinearForm<Scalar, Tag, Parameters...> operator-(
    SymmetricBilinearForm<Scalar, Tag, Parameters...> const& right);

template<quantity Scalar, typename Tag, typename... Parameters>
SymmetricBilinearForm<Scalar, Tag, Parameters...> operator+(
    SymmetricBilinearForm<Scalar, Tag, Parameters...> const& left,
    SymmetricBilinearForm<Scalar, Tag, Parameters...> const& right);
template<quantity Scalar, typename Tag, typename... Parameters>
SymmetricBilinearForm<Scalar, Tag, Parameters...> operator-(
    SymmetricBilinearForm<Scalar, Tag, Parameters...> const& left,
    SymmetricBilinearForm<Scalar, Tag, Parameters...> const& right);

// Vector space.

template<convertible_to_quantity LScalar,
         quantity RScalar,
         typename Tag,
         typename... Parameters>
SymmetricBilinearForm<Product<LScalar, RScalar>, Tag, Parameters...>
operator*(LScalar left,
          SymmetricBilinearForm<RScalar, Tag, Parameters...> const& right);
template<quantity LScalar,
         convertible_to_quantity RScalar,
         typename Tag,
         typename... Parameters>
SymmetricBilinearForm<Product<LScalar, RScalar>, Tag, Parameters...>
operator*(SymmetricBilinearForm<LScalar, Tag, Parameters...> const& left,
          RScalar right);
template<quantity LScalar,
         convertible_to_quantity RScalar,
         typename Tag,
         typename... Parameters>
SymmetricBilinearForm<Quotient<LScalar, RScalar>, Tag, Parameters...>
operator/(SymmetricBilinearForm<LScalar, Tag, Parameters...> const& left,
          RScalar right);

// For multivectors, it is convenient to write the application of a form as a
// product, because that's what is customary in physics.
// This only works because we identify vectors with linear forms on vectors; the
// application to a vector should return a linear form.  We also make that
// identification for quantities and reals (where it is canonical), but because
// we need to distinguish a symmetric bilinear form on the reals from a real,
// the `operator*` would be ambiguous.
template<quantity LScalar, quantity RScalar, typename Frame, int rank>
Multivector<Product<LScalar, RScalar>, Frame, rank> operator*(
    SymmetricBilinearForm<
        LScalar, MultivectorTag, Frame, Rank<rank>> const& left,
    Multivector<RScalar, Frame, rank> const& right);

template<quantity LScalar, quantity RScalar, typename Frame, int rank>
Multivector<Product<LScalar, RScalar>, Frame, rank> operator*(
    Multivector<LScalar, Frame, rank> const& left,
    SymmetricBilinearForm<
        RScalar, MultivectorTag, Frame, Rank<rank>> const& right);

// Solves the system `result * right == left`.  Note that by symmetry, this is
// also the solution of `right * result == left`.
template<quantity LScalar, quantity RScalar, typename Frame, int rank>
Multivector<Quotient<LScalar, RScalar>, Frame, rank> operator/(
    Multivector<LScalar, Frame, rank> const& left,
    SymmetricBilinearForm<
        RScalar, MultivectorTag, Frame, Rank<rank>> const& right);

// `SymmetricProduct(v, w)` is v ⊙ w ≔ (v ⊗ w + w ⊗ v) / 2.
template<quantity LScalar, quantity RScalar, typename Frame>
SymmetricBilinearForm<Product<LScalar, RScalar>, MultivectorTag, Frame, Rank<1>>
SymmetricProduct(Vector<LScalar, Frame> const& left,
                 Vector<RScalar, Frame> const& right);

template<quantity LScalar, quantity RScalar>
SymmetricBilinearForm<Product<LScalar, RScalar>, QuantityTag>
SymmetricProduct(LScalar const& left, RScalar const& right);

// `SymmetricSquare(v)` is `SymmetricProduct(v, v)`.
template<quantity Scalar, typename Frame>
SymmetricBilinearForm<Square<Scalar>, MultivectorTag, Frame, Rank<1>>
SymmetricSquare(Vector<Scalar, Frame> const& vector);

template<quantity Scalar>
SymmetricBilinearForm<Square<Scalar>, QuantityTag>
SymmetricSquare(Scalar const& quantity);

// Symmetric bilinear forms on vectors act on bivectors through this function.
// `Anticommutator(F, B)` is (tr(F)𝟙 - F)B in ℝ³ representation.  In matrix
// representation it is FB + BF = {F, B}.
template<quantity LScalar, quantity RScalar, typename Frame>
Bivector<Product<LScalar, RScalar>, Frame> Anticommutator(
    SymmetricBilinearForm<LScalar,
                          MultivectorTag, Frame, Rank<1>> const& form,
    Bivector<RScalar, Frame> const& bivector);

template<quantity Scalar, typename Tag, typename... Parameters>
std::string DebugString(
    SymmetricBilinearForm<Scalar, Tag, Parameters...> const& form);

template<quantity Scalar, typename Tag, typename... Parameters>
std::ostream& operator<<(
    std::ostream& out,
    SymmetricBilinearForm<Scalar, Tag, Parameters...> const& form);


// Convenient aliases.  Prefer using these aliases instead of the tag types.
template<typename Scalar, typename Frame>
using VectorSymmetricBilinearForm =
    SymmetricBilinearForm<Scalar, MultivectorTag, Frame, Rank<1>>;
template<typename Scalar, typename Frame>
using BivectorSymmetricBilinearForm =
    SymmetricBilinearForm<Scalar, MultivectorTag, Frame, Rank<2>>;
template<typename Scalar>
using QuantitySymmetricBilinearForm =
    SymmetricBilinearForm<Scalar, QuantityTag>;

}  // namespace internal

using internal::Anticommutator;
using internal::BivectorSymmetricBilinearForm;
using internal::QuantitySymmetricBilinearForm;
using internal::SymmetricBilinearForm;
using internal::SymmetricProduct;
using internal::SymmetricSquare;
using internal::VectorSymmetricBilinearForm;

}  // namespace _symmetric_bilinear_form
}  // namespace geometry
}  // namespace principia

#include "geometry/symmetric_bilinear_form_body.hpp"
