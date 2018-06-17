
#pragma once

#include "geometry/grassmann.hpp"
#include "quantities/named_quantities.hpp"

namespace principia {
namespace geometry {
namespace internal_symmetric_bilinear_form {

using quantities::Product;

template<typename Scalar, typename Frame>
class SymmetricBilinearForm {
 public:
  SymmetricBilinearForm();  //???

  //???dimensionful vector space structure

  template<typename LScalar, typename RScalar>
  Product<Scalar, Product<LScalar, RScalar>> operator()(
      Vector<LScalar, Frame> const& left,
      Vector<RScalar, Frame> const& right);

  template<typename LScalar, typename RScalar>
  Product<Scalar, Product<LScalar, RScalar>> operator()(
      Bivector<LScalar, Frame> const& left,
      Bivector<RScalar, Frame> const& right);

  static SymmetricBilinearForm InnerProductForm();
};

template<typename LScalar, typename RScalar, typename Frame>
Product<LScalar, RScalar> operator*(
    SymmetricBilinearForm<LScalar, Frame> const& left,
    Vector<RScalar, Frame> const& right);

template<typename LScalar, typename RScalar, typename Frame>
Product<LScalar, RScalar> operator*(
    SymmetricBilinearForm<LScalar, Frame> const& left,
    Bivector<RScalar, Frame> const& right);

template<typename LScalar, typename RScalar, typename Frame>
Product<LScalar, RScalar> operator*(
    Vector<LScalar, Frame> const& left,
    SymmetricBilinearForm<RScalar, Frame> const& right);

template<typename LScalar, typename RScalar, typename Frame>
Product<LScalar, RScalar> operator*(
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

}  // namespace internal_symmetric_bilinear_form

using internal_symmetric_bilinear_form::SymmetricBilinearForm;

}  // namespace geometry
}  // namespace principia

#include "geometry/symmetric_bilinear_form_body.hpp"
