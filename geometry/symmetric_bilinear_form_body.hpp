
#pragma once

#include "geometry/symmetric_bilinear_form.hpp"

namespace principia {
namespace geometry {
namespace internal_symmetric_bilinear_form {

template<typename Scalar, typename Frame>
SymmetricBilinearForm<Scalar, Frame>::SymmetricBilinearForm() {}

template<typename Scalar, typename Frame>
template<typename LScalar, typename RScalar>
Product<Scalar, Product<LScalar, RScalar>>
SymmetricBilinearForm<Scalar, Frame>::operator()(
    Vector<LScalar, Frame> const& left,
    Vector<RScalar, Frame> const& right) {
}

template<typename Scalar, typename Frame>
template<typename LScalar, typename RScalar>
Product<Scalar, Product<LScalar, RScalar>>
SymmetricBilinearForm<Scalar, Frame>::operator()(
    Bivector<LScalar, Frame> const& left,
    Bivector<RScalar, Frame> const& right) {
}

template<typename Scalar, typename Frame>
SymmetricBilinearForm<Scalar, Frame>
SymmetricBilinearForm<Scalar, Frame>::InnerProductForm() {}

template<typename LScalar, typename RScalar, typename Frame>
Product<LScalar, RScalar> operator*(
    SymmetricBilinearForm<LScalar, Frame> const& left,
    Vector<RScalar, Frame> const& right) {
}

template<typename LScalar, typename RScalar, typename Frame>
Product<LScalar, RScalar> operator*(
    SymmetricBilinearForm<LScalar, Frame> const& left,
    Bivector<RScalar, Frame> const& right) {
}

template<typename LScalar, typename RScalar, typename Frame>
Product<LScalar, RScalar> operator*(
    Vector<LScalar, Frame> const& left,
    SymmetricBilinearForm<RScalar, Frame> const& right) {
}

template<typename LScalar, typename RScalar, typename Frame>
Product<LScalar, RScalar> operator*(
    Bivector<LScalar, Frame> const& left,
    SymmetricBilinearForm<RScalar, Frame> const& right) {
}

template<typename LScalar, typename RScalar, typename Frame>
SymmetricBilinearForm<Product<LScalar, RScalar>, Frame> SymmetricProduct(
    Vector<LScalar, Frame> const& left,
    Vector<RScalar, Frame> const& right) {
}

template<typename LScalar, typename RScalar, typename Frame>
SymmetricBilinearForm<Product<LScalar, RScalar>, Frame> SymmetricProduct(
    Bivector<LScalar, Frame> const& left,
    Bivector<RScalar, Frame> const& right) {
}

}  // namespace internal_symmetric_bilinear_form
}  // namespace geometry
}  // namespace principia