#pragma once

#include "geometry/grassmann.hpp"
#include "geometry/symmetric_bilinear_form.hpp"
#include "quantities/named_quantities.hpp"

namespace principia {
namespace physics {
namespace _inertia_tensor {
namespace internal {

using namespace principia::geometry::_grassmann;
using namespace principia::geometry::_symmetric_bilinear_form;
using namespace principia::quantities::_named_quantities;

template<typename Frame>
using InertiaTensor = SymmetricBilinearForm<MomentOfInertia, Frame, Bivector>;

}  // namespace internal

using internal::InertiaTensor;

}  // namespace _inertia_tensor
}  // namespace physics
}  // namespace principia
