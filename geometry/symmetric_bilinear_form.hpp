
#pragma once

#include <string>

#include "base/macros.hpp"
#include "base/not_null.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/r3x3_matrix.hpp"
#include "quantities/named_quantities.hpp"
#include "serialization/geometry.pb.h"

namespace principia {
namespace geometry {

FORWARD_DECLARE_FROM(rotation,
                     TEMPLATE(typename FromFrame, typename ToFrame) class,
                     Rotation);

namespace internal_symmetric_bilinear_form {

using base::not_null;
using quantities::Product;
using quantities::Quotient;

template<typename Scalar, typename Frame>
class SymmetricBilinearForm {
 public:
  explicit SymmetricBilinearForm(R3x3Matrix<Scalar> const& matrix);
  explicit SymmetricBilinearForm(R3x3Matrix<Scalar>&& matrix);

  SymmetricBilinearForm& operator+=(SymmetricBilinearForm const& right);
  SymmetricBilinearForm& operator-=(SymmetricBilinearForm const& right);
  SymmetricBilinearForm& operator*=(double right);
  SymmetricBilinearForm& operator/=(double right);

  R3x3Matrix<Scalar> const& coordinates() const;

  // A SymmetricBilinearForm does *not* apply to bivectors.  See Anticommutator
  // below.
  template<typename LScalar, typename RScalar>
  Product<Scalar, Product<LScalar, RScalar>> operator()(
      Vector<LScalar, Frame> const& left,
      Vector<RScalar, Frame> const& right) const;

  // The eigensystem for a form is described by (1) the form in its eigenbasis,
  // which gives the eigenvalues; and (2) a rotation from the current basis to
  // the eigenbasis, which gives the eigenvectors.
  template<typename Eigenframe>
  struct Eigensystem {
    SymmetricBilinearForm<Scalar, Eigenframe> form;
    Rotation<Frame, Eigenframe> rotation;
  };

  // Computes a form equivalent to the current one but diagonalized with
  // increasing eigenvalues.
  template<typename Eigenframe>
  Eigensystem<Eigenframe> Diagonalize() const;

  void WriteToMessage(
      not_null<serialization::SymmetricBilinearForm*> message) const;
  template<typename = std::enable_if_t<base::is_serializable_v<Frame>>>
  static SymmetricBilinearForm ReadFromMessage(
      serialization::SymmetricBilinearForm const& message);

 private:
  // Given a matrix that contains in columns eigenvectors for a form, picks the
  // column with the largest norm and return its normalized value.  This is
  // useful to extract eigenvectors when eigenvalues are known.
  template<typename S>
  static R3Element<double> PickEigenvector(R3x3Matrix<S> const& matrix);

  // All the operations on this class must ensure that this matrix remains
  // symmetric.
  R3x3Matrix<Scalar> matrix_;

  template<typename S, typename F>
  friend class SymmetricBilinearForm;

  template<typename F>
  friend SymmetricBilinearForm<double, F> const& InnerProductForm();

  template<typename S, typename F>
  friend SymmetricBilinearForm<S, F> operator+(
      SymmetricBilinearForm<S, F> const& right);
  template<typename S, typename F>
  friend SymmetricBilinearForm<S, F> operator-(
      SymmetricBilinearForm<S, F> const& right);

  template<typename S, typename F>
  friend SymmetricBilinearForm<S, F> operator+(
      SymmetricBilinearForm<S, F> const& left,
      SymmetricBilinearForm<S, F> const& right);
  template<typename S, typename F>
  friend SymmetricBilinearForm<S, F> operator-(
      SymmetricBilinearForm<S, F> const& left,
      SymmetricBilinearForm<S, F> const& right);

  template<typename L, typename R, typename F>
  friend SymmetricBilinearForm<Product<L, R>, F> operator*(
      L left,
      SymmetricBilinearForm<R, F> const& right);
  template<typename L, typename R, typename F>
  friend SymmetricBilinearForm<Product<L, R>, F> operator*(
      SymmetricBilinearForm<L, F> const& left,
      R right);
  template<typename L, typename R, typename F>
  friend SymmetricBilinearForm<Quotient<L, R>, F> operator/(
      SymmetricBilinearForm<L, F> const& left,
      R right);

  template<typename L, typename R, typename F>
  friend Vector<Product<L, R>, F> operator*(
      SymmetricBilinearForm<L, F> const& left,
      Vector<R, F> const& right);

  template<typename L, typename R, typename F>
  friend Vector<Product<L, R>, F> operator*(
      Vector<L, F> const& left,
      SymmetricBilinearForm<R, F> const& right);

  template<typename L, typename R, typename F>
  friend SymmetricBilinearForm<Product<L, R>, F>
  SymmetricProduct(Vector<L, F> const& left,
                   Vector<R, F> const& right);

  template<typename L, typename R, typename F>
  friend Bivector<Product<L, R>, F> Anticommutator(
      SymmetricBilinearForm<L, F> const& form,
      Bivector<R, F> const& bivector);

  template<typename S, typename F>
  friend bool operator==(SymmetricBilinearForm<S, F> const& left,
                         SymmetricBilinearForm<S, F> const& right);
  template<typename S, typename F>
  friend bool operator!=(SymmetricBilinearForm<S, F> const& left,
                         SymmetricBilinearForm<S, F> const& right);

  template<typename S, typename F>
  friend std::string DebugString(SymmetricBilinearForm<S, F> const& form);

  template<typename S, typename F>
  friend std::ostream& operator<<(std::ostream& out,
                                  SymmetricBilinearForm<S, F> const& form);

  friend class SymmetricBilinearFormTest;
};

// |InnerProductForm()| is the symmetric bilinear form such that for all v and
// w, |InnerProductForm()(v, w) == InnerProduct(v, w)|.
template<typename Frame>
SymmetricBilinearForm<double, Frame> const& InnerProductForm();

template<typename Scalar, typename Frame>
SymmetricBilinearForm<Scalar, Frame> operator+(
    SymmetricBilinearForm<Scalar, Frame> const& right);
template<typename Scalar, typename Frame>
SymmetricBilinearForm<Scalar, Frame> operator-(
    SymmetricBilinearForm<Scalar, Frame> const& right);

template<typename Scalar, typename Frame>
SymmetricBilinearForm<Scalar, Frame> operator+(
    SymmetricBilinearForm<Scalar, Frame> const& left,
    SymmetricBilinearForm<Scalar, Frame> const& right);
template<typename Scalar, typename Frame>
SymmetricBilinearForm<Scalar, Frame> operator-(
    SymmetricBilinearForm<Scalar, Frame> const& left,
    SymmetricBilinearForm<Scalar, Frame> const& right);

template<typename LScalar, typename RScalar, typename Frame>
SymmetricBilinearForm<Product<LScalar, RScalar>, Frame> operator*(
    LScalar left,
    SymmetricBilinearForm<RScalar, Frame> const& right);
template<typename LScalar, typename RScalar, typename Frame>
SymmetricBilinearForm<Product<LScalar, RScalar>, Frame> operator*(
    SymmetricBilinearForm<LScalar, Frame> const& left,
    RScalar right);
template<typename LScalar, typename RScalar, typename Frame>
SymmetricBilinearForm<Quotient<LScalar, RScalar>, Frame> operator/(
    SymmetricBilinearForm<LScalar, Frame> const& left,
    RScalar right);

template<typename LScalar, typename RScalar, typename Frame>
Vector<Product<LScalar, RScalar>, Frame> operator*(
    SymmetricBilinearForm<LScalar, Frame> const& left,
    Vector<RScalar, Frame> const& right);

template<typename LScalar, typename RScalar, typename Frame>
Vector<Product<LScalar, RScalar>, Frame> operator*(
    Vector<LScalar, Frame> const& left,
    SymmetricBilinearForm<RScalar, Frame> const& right);

// |SymmetricProduct(v, w)| is v ⊙ w ≔ (v ⊗ w + w ⊗ v) / 2.
template<typename LScalar, typename RScalar, typename Frame>
SymmetricBilinearForm<Product<LScalar, RScalar>, Frame> SymmetricProduct(
    Vector<LScalar, Frame> const& left,
    Vector<RScalar, Frame> const& right);

// Symmetric bilinear forms act on bivectors through this function.
// |Anticommutator(F, B)| is (tr(F)𝟙 - F)B in ℝ³ representation .  In matrix
// representation it is FB + BF = {F, B}.
template<typename LScalar, typename RScalar, typename Frame>
Bivector<Product<LScalar, RScalar>, Frame> Anticommutator(
    SymmetricBilinearForm<LScalar, Frame> const& form,
    Bivector<RScalar, Frame> const& bivector);

template<typename Scalar, typename Frame>
bool operator==(SymmetricBilinearForm<Scalar, Frame> const& left,
                SymmetricBilinearForm<Scalar, Frame> const& right);
template<typename Scalar, typename Frame>
bool operator!=(SymmetricBilinearForm<Scalar, Frame> const& left,
                SymmetricBilinearForm<Scalar, Frame> const& right);

template<typename Scalar, typename Frame>
std::string DebugString(SymmetricBilinearForm<Scalar, Frame> const& form);

template<typename Scalar, typename Frame>
std::ostream& operator<<(std::ostream& out,
                         SymmetricBilinearForm<Scalar, Frame> const& form);

}  // namespace internal_symmetric_bilinear_form

using internal_symmetric_bilinear_form::InnerProductForm;
using internal_symmetric_bilinear_form::SymmetricBilinearForm;
using internal_symmetric_bilinear_form::SymmetricProduct;

}  // namespace geometry
}  // namespace principia

#include "geometry/symmetric_bilinear_form_body.hpp"
