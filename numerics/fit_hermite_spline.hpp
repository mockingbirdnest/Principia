#pragma once

#include <list>

#include "numerics/hermite3.hpp"

namespace principia {
namespace numerics {
namespace internal_fit_hermite_spline {

// Given |samples| for which the arguments, values, and derivatives can be
// obtained via the given functors, returns iterators for which an interpolating
// Hermite spline "barely" satisfies the tolerance. The result contains only
// right endpoints of interpolated intervals (i.e., the first interpolating
// polynomial is [begin, result.front()]).  It may not contain |--end|, if the
// resulting polynomial would have an error too small before the tolerance: in
// particular the result may be empty.  |Samples| must be have random access
// iterators.
template<typename Samples,
         typename GetArgument,
         typename GetValue,
         typename GetDerivative,
         typename ErrorType>
std::list<typename Samples::iterator> FitHermiteSpline(
    Samples const& samples,
    GetArgument get_argument,
    GetValue get_value,
    GetDerivative get_derivative,
    ErrorType tolerance);

}  // namespace internal_fit_hermite_spline

using internal_fit_hermite_spline::FitHermiteSpline;

}  // namespace numerics
}  // namespace principia

#include "numerics/fit_hermite_spline_body.hpp"
