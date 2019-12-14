#pragma once

#include "physics/inertia_tensor.hpp"

#include "geometry/barycentre_calculator.hpp"
#include "geometry/identity.hpp"
#include "geometry/point.hpp"
#include "geometry/rotation.hpp"
#include "geometry/symmetric_bilinear_form.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace physics {
namespace internal_inertia_tensor {

using geometry::Anticommutator;
using geometry::Barycentre;
using geometry::Bivector;
using geometry::Commutator;
using geometry::Identity;
using geometry::Rotation;
using geometry::SymmetricProduct;
using quantities::Cbrt;
using quantities::Density;
using quantities::Pow;
using quantities::si::Kilogram;
using quantities::si::Metre;

namespace {
constexpr MomentOfInertia zero;
}  // namespace

template<typename Frame>
InertiaTensor<Frame>::InertiaTensor()
    : InertiaTensor(Mass{},
                    R3x3Matrix<MomentOfInertia>{{zero, zero, zero},
                                                {zero, zero, zero},
                                                {zero, zero, zero}},
                    Frame::origin) {}

template<typename Frame>
InertiaTensor<Frame>::InertiaTensor(
    Mass const& mass,
    R3x3Matrix<MomentOfInertia> const& coordinates,
    Position<Frame> const& centre_of_mass)
    : InertiaTensor(mass,
                    MakeSymmetricBilinearForm(coordinates),
                    centre_of_mass) {}

template<typename Frame>
Mass const& InertiaTensor<Frame>::mass() const {
  return mass_;
}

template<typename Frame>
Position<Frame> const& InertiaTensor<Frame>::centre_of_mass() const {
  return centre_of_mass_;
}

template<typename Frame>
template<typename ToFrame>
InertiaTensor<ToFrame> InertiaTensor<Frame>::Transform(
    RigidTransformation<Frame, ToFrame> const& transformation) const {
  Displacement<ToFrame> const displacement =
      ToFrame::origin - transformation(Frame::origin);
  SymmetricBilinearForm<MomentOfInertia, ToFrame> const transformed_form =
      transformation.linear_map()(form_) +
      mass_ * SymmetricProduct(displacement, displacement);
  Position<ToFrame> const transformed_centre_of_mass =
      transformation(centre_of_mass_);
  return InertiaTensor<ToFrame>(mass_,
                                transformed_form,
                                transformed_centre_of_mass);
}

template<typename Frame>
template<typename PrincipalAxesFrame>
typename InertiaTensor<Frame>::template PrincipalAxes<PrincipalAxesFrame>
InertiaTensor<Frame>::Diagonalize() const {
  struct IntermediateFrame {};

  auto const eigensystem = form_.template Diagonalize<IntermediateFrame>();

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

  auto const y = Bivector<double, IntermediateFrame>({0, 1, 0});
  auto const z = Bivector<double, IntermediateFrame>({0, 0, 1});
  Rotation<IntermediateFrame, PrincipalAxesFrame> const swap_x_and_z(
      z, y, Commutator(z, y));

  return {moments_of_inertia, swap_x_and_z * eigensystem.rotation};
}

template<typename Frame>
InertiaTensor<Frame> InertiaTensor<Frame>::MakeWaterSphereInertiaTensor(
    Mass const& mass) {
  static constexpr Density ρ_of_water = 1000 * Kilogram / Pow<3>(Metre);
  MomentOfInertia const I =
      Cbrt(9 * Pow<5>(mass) / (250 * Pow<2>(π * ρ_of_water)));
  return InertiaTensor<Frame>(mass,
                              R3x3Matrix<MomentOfInertia>({I, zero, zero},
                                                          {zero, I, zero},
                                                          {zero, zero, I}),
                              Frame::origin);
}

template<typename Frame>
void InertiaTensor<Frame>::WriteToMessage(
    not_null<serialization::InertiaTensor*> const message) const {
  mass_.WriteToMessage(message->mutable_mass());
  form_.WriteToMessage(message->mutable_form());
  centre_of_mass_.WriteToMessage(message->mutable_centre_of_mass());
}

template<typename Frame>
InertiaTensor<Frame> InertiaTensor<Frame>::ReadFromMessage(
    serialization::InertiaTensor const& message) {
  return InertiaTensor(
      Mass::ReadFromMessage(message.mass()),
      SymmetricBilinearForm<MomentOfInertia, Frame>::ReadFromMessage(
          message.form()),
      Position<Frame>::ReadFromMessage(message.centre_of_mass()));
}

template<typename Frame>
InertiaTensor<Frame>::InertiaTensor(
    Mass const& mass,
    SymmetricBilinearForm<MomentOfInertia, Frame> const& form,
    Position<Frame> const& centre_of_mass)
    : mass_(mass),
      form_(form),
      centre_of_mass_(centre_of_mass) {
  CHECK_GE(mass, Mass{});
}

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
           0.5 * (tensor(0, 0) + tensor(1, 1) - tensor(2, 2))}));
}

template<typename RScalar, typename Frame>
Bivector<Product<MomentOfInertia, RScalar>, Frame> Anticommutator(
    InertiaTensor<Frame> const& tensor,
    Bivector<RScalar, Frame> const& bivector) {
  return Anticommutator(tensor.form_, bivector);
}

template<typename Frame>
bool operator==(InertiaTensor<Frame> const& left,
                InertiaTensor<Frame> const& right) {
  return left.mass_ == right.mass_ &&
         left.form_ == right.form_ &&
         left.centre_of_mass_ == right.centre_of_mass_;
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

template<typename Frame>
InertiaTensor<Frame>& operator+=(InertiaTensor<Frame>& left,
                                 InertiaTensor<Frame> const& right) {
  left = left + right;
  return left;
}

}  // namespace internal_inertia_tensor
}  // namespace physics
}  // namespace principia
