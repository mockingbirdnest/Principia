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
// It must declare the type of the gradient and the inverse Hessian thus:
//   using Gradient = ...;
//   using InverseHessian = ...;
template<homogeneous_field Scalar, homogeneous_euclidean_space<Scalar> Argument>
struct Generator;

// The type of the gradient of `VectorField<Scalar, Argument>`.
template<homogeneous_field Scalar, homogeneous_euclidean_space<Scalar> Argument>
using Gradient = typename Generator<Scalar, Argument>::Gradient;

// The gradient vector field.
template<homogeneous_field Scalar, homogeneous_euclidean_space<Scalar> Argument>
using GradientField = VectorField<Gradient<Scalar, Argument>, Argument>;

// The type of the inverse Hessian of `VectorField<Scalar, Argument>`.
template<homogeneous_field Scalar, homogeneous_euclidean_space<Scalar> Argument>
using InverseHessian = typename Generator<Scalar, Argument>::InverseHessian;

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
using internal::InverseHessian;
using internal::VectorField;

}  // namespace _differential_geometry
}  // namespace geometry
}  // namespace principia

#include "geometry/differential_geometry_body.hpp"
