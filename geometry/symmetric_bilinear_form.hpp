
#pragma once

#include "geometry/grassmann.hpp"
#include "geometry/r3x3_matrix.hpp"
#include "quantities/named_quantities.hpp"

namespace principia {
namespace geometry {
namespace internal_symmetric_bilinear_form {

using quantities::Product;

template<typename Scalar, typename Frame>
class SymmetricBilinearForm {
 public:
  SymmetricBilinearForm& operator+=(SymmetricBilinearForm const& right);
  SymmetricBilinearForm& operator-=(SymmetricBilinearForm const& right);
  SymmetricBilinearForm& operator*=(double right);
  SymmetricBilinearForm& operator/=(double right);

  template<typename LScalar, typename RScalar>
  Product<Scalar, Product<LScalar, RScalar>> operator()(
      Vector<LScalar, Frame> const& left,
      Vector<RScalar, Frame> const& right);

  template<typename LScalar, typename RScalar>
  Product<Scalar, Product<LScalar, RScalar>> operator()(
      Bivector<LScalar, Frame> const& left,
      Bivector<RScalar, Frame> const& right);

  static SymmetricBilinearForm InnerProductForm();

 private:
  explicit SymmetricBilinearForm(R3x3Matrix<Scalar> const& matrix);
  explicit SymmetricBilinearForm(R3x3Matrix<Scalar>&& matrix);

  // All the operations on this class must ensure that this matrix remains
  // symmetric.
  R3x3Matrix<Scalar> matrix_;
};

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

template<typename Scalar, typename Frame>
SymmetricBilinearForm<Scalar, Frame> operator*(
    double left,
    SymmetricBilinearForm<Scalar, Frame> const& right);
template<typename Scalar, typename Frame>
SymmetricBilinearForm<Scalar, Frame> operator*(
    SymmetricBilinearForm<Scalar, Frame> const& left,
    double right);
template<typename Scalar, typename Frame>
SymmetricBilinearForm<Scalar, Frame> operator/(
    SymmetricBilinearForm<Scalar, Frame> const& left,
    double right);

template<typename LScalar, typename RScalar, typename Frame>
Vector<Product<LScalar, RScalar>, Frame> operator*(
    SymmetricBilinearForm<LScalar, Frame> const& left,
    Vector<RScalar, Frame> const& right);

template<typename LScalar, typename RScalar, typename Frame>
Bivector<Product<LScalar, RScalar>, Frame> operator*(
    SymmetricBilinearForm<LScalar, Frame> const& left,
    Bivector<RScalar, Frame> const& right);

template<typename LScalar, typename RScalar, typename Frame>
Vector<Product<LScalar, RScalar>, Frame> operator*(
    Vector<LScalar, Frame> const& left,
    SymmetricBilinearForm<RScalar, Frame> const& right);

template<typename LScalar, typename RScalar, typename Frame>
Bivector<Product<LScalar, RScalar>, Frame> operator*(
    Bivector<LScalar, Frame> const& left,
    SymmetricBilinearForm<RScalar, Frame> const& right);

template<typename LScalar, typename RScalar, typename Frame>
SymmetricBilinearForm<Product<LScalar, RScalar>, Frame> SymmetricProduct(
    Vector<LScalar, Frame> const& left,
    Vector<RScalar, Frame> const& right);

template<typename LScalar, typename RScalar, typename Frame>
SymmetricBilinearForm<Product<LScalar, RScalar>, Frame> SymmetricProduct(
    Bivector<LScalar, Frame> const& left,
    Bivector<RScalar, Frame> const& right);

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

using internal_symmetric_bilinear_form::SymmetricBilinearForm;
using internal_symmetric_bilinear_form::SymmetricProduct;

}  // namespace geometry
}  // namespace principia

#include "geometry/symmetric_bilinear_form_body.hpp"
