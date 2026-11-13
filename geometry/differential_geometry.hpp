#pragma once

#include "base/algebra.hpp"
#include "geometry/hilbert.hpp"

namespace principia {
namespace geometry {
namespace _differential_geometry {
namespace internal {

using namespace principia::base::_algebra;
using namespace principia::geometry::_hilbert;

template<real_vector_space Value, euclidean Argument>
using VectorField = std::function<Value(Argument const&)>;

// A helper for generating declarations derived from `Scalar` and `Argument`.
// It must declare the type of the gradient and a function returning the inner
// product form thus:
//   using Gradient = ...;
//   static ... InnerProductForm();
// TODO(phl): The `InnerProductForm` should probably go elsewhere.
template<homogeneous_field Scalar, homogeneous_euclidean_space<Scalar> Argument>
struct Generator;

// The type of the gradient of `VectorField<Scalar, Argument>`.
template<homogeneous_field Scalar, homogeneous_euclidean_space<Scalar> Argument>
using Gradient = typename Generator<Scalar, Argument>::Gradient;

// The gradient vector field.
template<homogeneous_field Scalar, homogeneous_euclidean_space<Scalar> Argument>
using GradientField = VectorField<Gradient<Scalar, Argument>, Argument>;

// Must return InnerProduct(grad_f(argument), direction).  This is not a
// `VectorField` because it depends on the direction.
template<homogeneous_field Scalar, homogeneous_euclidean_space<Scalar> Argument>
using GateauxDerivative =
    std::function<Scalar(Argument const& argument,
                         Difference<Argument> const& direction)>;

}  // namespace internal

using internal::GateauxDerivative;
using internal::Gradient;
using internal::GradientField;
using internal::VectorField;

}  // namespace _differential_geometry
}  // namespace geometry
}  // namespace principia

#include "geometry/differential_geometry_body.hpp"
