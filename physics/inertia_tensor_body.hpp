#pragma once

#include "physics/inertia_tensor.hpp"

#include "geometry/barycentre_calculator.hpp"
#include "geometry/identity.hpp"
#include "geometry/point.hpp"
#include "geometry/rotation.hpp"

namespace principia {
namespace physics {
namespace internal_inertia_tensor {

using geometry::Barycentre;
using geometry::Bivector;
using geometry::Identity;
using geometry::Rotation;
using geometry::SymmetricProduct;

template<typename Frame>
InertiaTensor<Frame>::InertiaTensor(
    Mass const& mass,
    R3x3Matrix<MomentOfInertia> const& coordinates,
    Position<Frame> const& centre_of_mass)
    : InertiaTensor(mass,
                    MakeSymmetricBilinearForm(coordinates),
                    centre_of_mass) {}

template<typename Frame>
Position<Frame> const& InertiaTensor<Frame>::centre_of_mass() const {
  return centre_of_mass_;
}

template<typename Frame>
template<typename ToFrame>
InertiaTensor<ToFrame> InertiaTensor<Frame>::Rotate(
    Rotation<Frame, ToFrame> const& rotation) const {
  auto const centre_of_mass_displacement = centre_of_mass_ - Frame::origin;
  auto const rotated_centre_of_mass_displacement =
      rotation(centre_of_mass_displacement);
  return InertiaTensor<ToFrame>(
      mass_,
      rotation(form_),
      ToFrame::origin + rotated_centre_of_mass_displacement);
}

template<typename Frame>
template<typename ToFrame>
InertiaTensor<ToFrame> InertiaTensor<Frame>::Rotate(
    Rotation<ToFrame, Frame> const& rotation) const {
  return Rotate(rotation.Inverse());
}

template<typename Frame>
template<typename ToFrame>
InertiaTensor<ToFrame> InertiaTensor<Frame>::Translate(
    Position<Frame> const& point) const {
  static Identity<Frame, ToFrame> const identity{};

  auto const translation = point - Frame::origin;
  auto const translated_form =
      form_ + mass_ * SymmetricProduct(translation, translation);
  auto const translated_centre_of_mass_displacement = centre_of_mass_ - point;
  return InertiaTensor<ToFrame>(
      mass_,
      identity(translated_form),
      ToFrame::origin + identity(translated_centre_of_mass_displacement));
}

template<typename Frame>
template<typename PrincipalAxesFrame>
typename InertiaTensor<Frame>::PrincipalAxes<PrincipalAxesFrame>
InertiaTensor<Frame>::Diagonalize() const {
  struct IntermediateFrame {};

  auto const eigensystem = form_.Diagonalize<IntermediateFrame>();

  // The eigenvalues are {Σx², Σy², Σz²} in increasing order.
  R3x3Matrix<MomentOfInertia> const eigensystem_form_coordinates =
      eigensystem.form.coordinates();

  // The moment of inertia in increasing order are
  // {Σ(x² + y²), Σ(x² + z²), Σ(y² + z²)}, which does *not* correspond to the
  // ordered eigenvalues of the "classical" inertia tensor, which are
  // {Σ(y² + z²), Σ(x² + z²), Σ(x² + y²)}.  We address this by multiplying the
  // rotation by another one that swaps x and z.
  R3Element<MomentOfInertia> const moments_of_inertia(
      eigensystem_form_coordinates(0, 0) + eigensystem_form_coordinates(1, 1),
      eigensystem_form_coordinates(2, 2) + eigensystem_form_coordinates(0, 0),
      eigensystem_form_coordinates(1, 1) + eigensystem_form_coordinates(2, 2));

  Rotation<IntermediateFrame, PrincipalAxesFrame> const swap_x_and_z(
      Bivector<double, IntermediateFrame>({0, 0, 1}),
      Bivector<double, IntermediateFrame>({0, 1, 0}),
      Bivector<double, IntermediateFrame>({-1, 0, 0}));

  return {moments_of_inertia, swap_x_and_z * eigensystem.rotation};
}

template<typename Frame>
InertiaTensor<Frame>::InertiaTensor(
    Mass const& mass,
    SymmetricBilinearForm<MomentOfInertia, Frame> const& form,
    Position<Frame> const& centre_of_mass)
    : mass_(mass),
      form_(form),
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
  return InertiaTensor<Frame>(left.mass_ + right.mass_,
                              left.form_ + right.form_,
                              Barycentre<Position<Frame>, Mass>(
                                  {left.centre_of_mass_, right.centre_of_mass_},
                                  {left.mass_, right.mass_}));
}

}  // namespace internal_inertia_tensor
}  // namespace physics
}  // namespace principia