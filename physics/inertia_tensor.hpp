#pragma once

#include "base/not_null.hpp"
#include "geometry/named_quantities.hpp"
#include "geometry/point.hpp"
#include "geometry/r3_element.hpp"
#include "geometry/r3x3_matrix.hpp"
#include "geometry/rotation.hpp"
#include "geometry/symmetric_bilinear_form.hpp"
#include "physics/rigid_motion.hpp"
#include "quantities/named_quantities.hpp"
#include "serialization/physics.pb.h"

namespace principia {
namespace physics {
namespace internal_inertia_tensor {

using base::not_null;
using geometry::Displacement;
using geometry::Position;
using geometry::R3Element;
using geometry::R3x3Matrix;
using geometry::Rotation;
using geometry::SymmetricBilinearForm;
using quantities::Mass;
using quantities::MomentOfInertia;

// The inertia tensor of a solid at the origin and in the axis of Frame.
template<typename Frame>
class InertiaTensor {
 public:
  // Constructs an inertia tensor with the specified coordinates for an object
  // of the given mass.  The centre of mass of the object must be specified.
  // Typically, this constructor is used to implement formulæ like those in
  // https://en.wikipedia.org/wiki/List_of_moments_of_inertia and the centre of
  // mass is often Frame::origin.
  InertiaTensor(Mass const& mass,
                R3x3Matrix<MomentOfInertia> const& coordinates,
                Position<Frame> const& centre_of_mass);

  // The mass passed at construction.
  Mass const& mass() const;

  // Obtains the centre of mass of the object in Frame.
  Position<Frame> const& centre_of_mass() const;

  // Computes the inertia tensor in ToFrame, which is transformed from Frame
  // using the given rigid transformation.
  template<typename ToFrame>
  InertiaTensor<ToFrame> Transform(
      RigidTransformation<Frame, ToFrame> const& transformation) const;

  template<typename PrincipalAxesFrame>
  struct PrincipalAxes {
    R3Element<MomentOfInertia> moments_of_inertia;
    Rotation<Frame, PrincipalAxesFrame> rotation;
  };

  // A factory that creates an inertia tensor for a solid sphere of water having
  // the given mass.  Useful e.g. for save compatibility.
  static InertiaTensor<Frame> MakeWaterSphereInertiaTensor(Mass const& mass);

  // Diagonalization is possible in any frame, but it's mostly used in a frame
  // centred at the centre of mass.
  template<typename PrincipalAxesFrame>
  PrincipalAxes<PrincipalAxesFrame> Diagonalize() const;

  void WriteToMessage(not_null<serialization::InertiaTensor*> message) const;
  static InertiaTensor ReadFromMessage(
      serialization::InertiaTensor const& message);

 private:
  // Important: the form used internally to represent the inertia tensor does
  // *not* follow the convention customary in physics.  The usual convention is
  //  ⎛ ∑(y² + z²)   -∑xy      -∑xy    ⎞
  //  ⎜    -∑yx   ∑(x² + z²)   -∑yz    ⎟
  //  ⎝    -∑zx      -∑zy   ∑(x² + y²) ⎠
  // but our convention is:
  //  ⎛ ∑x²   ∑xy   ∑xy ⎞
  //  ⎜ ∑yx   ∑y²   ∑yz ⎟
  //  ⎝ ∑zx   ∑zy   ∑z² ⎠
  // This leads to simpler transformation formulæ at the expense of some care
  // when actually computing moments of inertia.
  InertiaTensor(Mass const& mass,
                SymmetricBilinearForm<MomentOfInertia, Frame> const& form,
                Position<Frame> const& centre_of_mass);

  // Computes our representation from the conventional representation.
  static SymmetricBilinearForm<MomentOfInertia, Frame>
  MakeSymmetricBilinearForm(R3x3Matrix<MomentOfInertia> const& tensor);

  Mass mass_;
  SymmetricBilinearForm<MomentOfInertia, Frame> form_;
  Position<Frame> centre_of_mass_;

  template<typename F>
  friend class InertiaTensor;

  template<typename F>
  friend bool operator==(InertiaTensor<F> const& left,
                         InertiaTensor<F> const& right);
  template<typename F>
  friend InertiaTensor<F> operator+(InertiaTensor<F> const& left,
                                    InertiaTensor<F> const& right);

  friend class InertiaTensorTest;
};

template<typename Frame>
bool operator==(InertiaTensor<Frame> const& left,
                InertiaTensor<Frame> const& right);

template<typename Frame>
InertiaTensor<Frame> operator+(InertiaTensor<Frame> const& left,
                               InertiaTensor<Frame> const& right);

}  // namespace internal_inertia_tensor

using internal_inertia_tensor::InertiaTensor;
using internal_inertia_tensor::operator==;
using internal_inertia_tensor::operator+;

}  // namespace physics
}  // namespace principia

#include "physics/inertia_tensor_body.hpp"
