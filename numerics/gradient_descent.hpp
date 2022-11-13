#pragma once

#include <functional>

#include "geometry/hilbert.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"

namespace principia {
namespace numerics {
namespace internal_gradient_descent {

using geometry::Hilbert;
using quantities::Difference;
using quantities::Length;
using quantities::Product;
using quantities::Quotient;

// In this file |Argument| must be such that its difference belongs to a Hilbert
// space.

template<typename Scalar, typename Argument>
using Field = std::function<Scalar(Argument const&)>;

template<typename Scalar, typename Argument>
using Gradient =
    Product<Scalar,
            Quotient<Difference<Argument>,
                     typename Hilbert<Difference<Argument>>::Norm²Type>>;

template<typename Scalar, typename Argument>
Argument BroydenFletcherGoldfarbShanno(
    Argument const& start_argument,
    Field<Scalar, Argument> const& f,
    Field<Gradient<Scalar, Argument>, Argument> const& grad_f,
    Length const& tolerance);

}  // namespace internal_gradient_descent

using internal_gradient_descent::BroydenFletcherGoldfarbShanno;

}  // namespace numerics
}  // namespace principia

#include "numerics/gradient_descent_body.hpp"
