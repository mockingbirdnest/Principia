#pragma once

#include "geometry/grassmann.hpp"
#include "geometry/symmetric_bilinear_form.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"

namespace principia {
namespace physics {
namespace _tensors {
namespace internal {

using namespace principia::geometry::_grassmann;
using namespace principia::geometry::_symmetric_bilinear_form;
using namespace principia::quantities::_named_quantities;
using namespace principia::quantities::_quantities;

template<typename Frame>
using InertiaTensor = SymmetricBilinearForm<MomentOfInertia, Frame, Bivector>;

// The type of the Jacobian of the gravitational acceleration field.
template<typename Frame>
using JacobianOfAcceleration =
    SymmetricBilinearForm<Inverse<Square<Time>>, Frame, Vector>;

}  // namespace internal

using internal::InertiaTensor;
using internal::JacobianOfAcceleration;

}  // namespace _tensors
}  // namespace physics
}  // namespace principia
