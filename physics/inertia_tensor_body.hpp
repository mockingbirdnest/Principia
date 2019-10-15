#pragma once

#include "physics/inertia_tensor.hpp"

namespace principia {
namespace physics {
namespace internal_inertia_tensor {

template<typename Frame>
InertiaTensor<Frame>::InertiaTensor(
    R3x3Matrix<MomentOfInertia> const& coordinates,
    Point<Frame> const& centre_of_mass) {}

template<typename Frame>
InertiaTensor<Frame>::InertiaTensor(
    R3x3Matrix<MomentOfInertia> const& coordinates,
    Point<Frame> const& reference_point,
    Point<Frame> const& centre_of_mass) {}

template<typename Frame>
R3Element<MomentOfInertia> InertiaTensor<Frame>::MomentsOfInertia() const {
  return R3Element<MomentOfInertia>();
}

template<typename Frame>
template<typename ToFrame>
InertiaTensor<ToFrame> InertiaTensor<Frame>::Rotate(
    Rotation<Frame, ToFrame> const& rotation) {
  return InertiaTensor<ToFrame>();
}

template<typename Frame>
template<typename ToFrame>
InertiaTensor<ToFrame> InertiaTensor<Frame>::Translate(
    Displacement<Frame> const& displacement) {
  return InertiaTensor<ToFrame>();
}

template<typename Frame>
template<typename ToFrame>
InertiaTensor<ToFrame> InertiaTensor<Frame>::Translate(
    Point<Frame> const& point) {
  return InertiaTensor<ToFrame>();
}

template<typename Frame>
template<typename PrincipalAxesFrame>
typename InertiaTensor<Frame>::PrincipalAxes<PrincipalAxesFrame>
InertiaTensor<Frame>::Diagonalize() const {
  return PrincipalAxes<PrincipalAxesFrame>();
}


}  // namespace internal_inertia_tensor
}  // namespace physics
}  // namespace principia