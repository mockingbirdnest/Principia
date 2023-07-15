#pragma once

#include <functional>
#include <optional>

#include "geometry/grassmann.hpp"
#include "geometry/hilbert.hpp"
#include "geometry/point.hpp"
#include "geometry/symmetric_bilinear_form.hpp"
#include "numerics/fixed_arrays.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"

namespace principia {
namespace numerics {
namespace _gradient_descent {
namespace internal {

using namespace principia::geometry::_grassmann;
using namespace principia::geometry::_hilbert;
using namespace principia::geometry::_point;
using namespace principia::geometry::_symmetric_bilinear_form;
using namespace principia::numerics::_fixed_arrays;
using namespace principia::quantities::_named_quantities;
using namespace principia::quantities::_quantities;

// In this file |Argument| must be such that its difference belongs to a Hilbert
// space.

// A helper for generating declarations derived from |Scalar| and |Argument|.
// It must declare the type of the gradient and a function returning the inner
// product form thus:
//   using Gradient = ...;
//   static ... InnerProductForm();
template<typename Scalar, typename Argument>
struct Generator;

template<typename Scalar, typename Argument>
using Field = std::function<Scalar(Argument const&)>;

template<typename Scalar, typename Argument>
using Gradient = typename Generator<Scalar, Argument>::Gradient;

// Stops when the search displacement is smaller than |tolerance|.  Returns
// |nullopt| if no minimum is found within distance |radius| of
// |start_argument|.
template<typename Scalar, typename Argument>
std::optional<Argument> BroydenFletcherGoldfarbShanno(
    Argument const& start_argument,
    Field<Scalar, Argument> const& f,
    Field<Gradient<Scalar, Argument>, Argument> const& grad_f,
    typename Hilbert<Difference<Argument>>::NormType const& tolerance,
    typename Hilbert<Difference<Argument>>::NormType const& radius =
        Infinity<typename Hilbert<Difference<Argument>>::NormType>);

}  // namespace internal

using internal::BroydenFletcherGoldfarbShanno;

}  // namespace _gradient_descent
}  // namespace numerics
}  // namespace principia

#include "numerics/gradient_descent_body.hpp"
