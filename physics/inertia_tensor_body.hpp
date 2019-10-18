#pragma once

#include "physics/inertia_tensor.hpp"

#include "geometry/barycentre_calculator.hpp"
#include "geometry/point.hpp"

namespace principia {
namespace physics {
namespace internal_inertia_tensor {

using geometry::Barycentre;
using geometry::SymmetricProduct;

template<typename Frame>
InertiaTensor<Frame>::InertiaTensor(
    Mass const& mass,
    R3x3Matrix<MomentOfInertia> const& coordinates,
    Position<Frame> const& centre_of_mass)
    : InertiaTensor(mass, coordinates, centre_of_mass, centre_of_mass) {}

template<typename Frame>
InertiaTensor<Frame>::InertiaTensor(
    Mass const& mass,
    R3x3Matrix<MomentOfInertia> const& coordinates,
    Position<Frame> const& reference_point,
    Position<Frame> const& centre_of_mass)
    : InertiaTensor(mass,
                    MakeSymmetricBilinearForm(coordinates),
                    reference_point,
                    centre_of_mass) {}

template<typename Frame>
R3Element<MomentOfInertia> InertiaTensor<Frame>::MomentsOfInertia() const {
  return R3Element<MomentOfInertia>();
}

template<typename Frame>
template<typename ToFrame>
InertiaTensor<ToFrame> InertiaTensor<Frame>::Rotate(
    Rotation<Frame, ToFrame> const& rotation) const {
  return InertiaTensor<ToFrame>();
}

template<typename Frame>
template<typename ToFrame>
InertiaTensor<ToFrame> InertiaTensor<Frame>::Translate(
    Displacement<Frame> const& displacement) const {
  return InertiaTensor(
      mass_,
      form_ + 2 * mass_ * SymmetricProduct(displacement, displacement),
      reference_point_ + displacement,
      centre_of_mass_);
}

template<typename Frame>
template<typename ToFrame>
InertiaTensor<ToFrame> InertiaTensor<Frame>::Translate(
    Displacement<ToFrame> const& displacement) const {
  return InertiaTensor(
      mass_,
      form_ + 2 * mass_ * SymmetricProduct(displacement, displacement),
      reference_point_ + displacement,
      centre_of_mass_);
}

template<typename Frame>
template<typename ToFrame>
InertiaTensor<ToFrame> InertiaTensor<Frame>::Translate(
    Position<Frame> const& point) const {
  auto const displacement = point - reference_point_;
  return InertiaTensor(
      mass_,
      form_ + 2 * mass_ * SymmetricProduct(displacement, displacement),
      point,
      centre_of_mass_);
}

template<typename Frame>
template<typename ToFrame>
InertiaTensor<ToFrame> InertiaTensor<Frame>::Translate(
    Position<ToFrame> const& point) const {
  auto const displacement = point - reference_point_;
  return InertiaTensor(
      mass_,
      form_ + 2 * mass_ * SymmetricProduct(displacement, displacement),
      point,
      centre_of_mass_);
}

template<typename Frame>
template<typename PrincipalAxesFrame>
typename InertiaTensor<Frame>::PrincipalAxes<PrincipalAxesFrame>
InertiaTensor<Frame>::Diagonalize() const {
  auto const eigensystem = form_.Diagonalize();
  //TODO(phl): What happens to the points?  reference_point_: identity;
  // centre_of_mass_: rotation.
  return {InertiaTensor<PrincipalAxesFrame>(
              mass_, eigensystem.form, reference_point_, centre_of_mass_),
          eigensystem.rotation};
}

template<typename Frame>
InertiaTensor<Frame>::InertiaTensor(
    Mass const& mass,
    SymmetricBilinearForm<MomentOfInertia, Frame> const& form,
    Position<Frame> const& reference_point,
    Position<Frame> const& centre_of_mass)
    : mass_(mass),
      form_(form),
      reference_point_(reference_point),
      centre_of_mass_(centre_of_mass) {}

template<typename Frame>
SymmetricBilinearForm<MomentOfInertia, Frame>
InertiaTensor<Frame>::MakeSymmetricBilinearForm(
    R3x3Matrix<MomentOfInertia> const& tensor) {
  return SymmetricBilinearForm<MomentOfInertia, Frame>(
      R3x3Matrix<MomentOfInertia>(
          {0.5 * (tensor(1, 1) + tensor(2, 2) - tensor(0, 0)),
           -tensor(0, 1),
           -tensor(0, 2)},
          {-tensor(1, 0),
           0.5 * (tensor(2, 2) + tensor(0, 0) - tensor(1, 1)),
           -tensor(1, 2)},
          {-tensor(2, 0),
           -tensor(2, 1),
           0.5 * (tensor(0, 0) + tensor(1, 1) - tensor(2, 2)) / 2}));
}

template<typename Frame>
InertiaTensor<Frame> operator+(InertiaTensor<Frame> const& left,
                               InertiaTensor<Frame> const& right) {
  CHECK_EQ(left.reference_point_, right.reference_point_);
  return InertiaTensor<Frame>(left.mass_ + right.mass_,
                              left.form_ + right.form_,
                              left.reference_point_,
                              Barycentre<Position<Frame>, Mass>(
                                  {left.centre_of_mass_, right.centre_of_mass_},
                                  {left.mass_, right.mass_}));
}

}  // namespace internal_inertia_tensor
}  // namespace physics
}  // namespace principia