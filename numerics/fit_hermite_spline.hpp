#pragma once

#include <list>

#include "geometry/hilbert.hpp"
#include "numerics/hermite3.hpp"

namespace principia {
namespace numerics {
namespace internal_fit_hermite_spline {

using quantities::Derivative;
using quantities::Difference;
using geometry::Hilbert;

// Given |samples| for which the arguments, values, and derivatives can be
// obtained via the given functors, returns a sequence of iterators
// (it₀, it₁, it₂, ..., itᵣ) such that, with it₋₁ = samples.begin(),
// for all integers it in [0, r], the |Hermite3| interpolation of
// (*itᵢ₋₁, *itᵢ) fits |samples| within |tolerance|, but the interpolation
// of (*itᵢ₋₁, *(it + 1)ᵢ) does not, i.e., the iterators delimit "maximal"
// fitting polynomials.
// Note that it follows from the above that itᵣ < samples.end() - 1, so that at
// not all of |samples| is fitted maximally.  This function further guarantees
// that the |Hermite3| interpolation of (*itᵣ, *(samples.end() - 1)) fits
// |samples| within |tolerance|.
template<typename Argument, typename Value, typename Samples>
std::list<typename Samples::const_iterator> FitHermiteSpline(
    Samples const& samples,
    std::function<Argument const&(typename Samples::value_type const&)> const&
        get_argument,
    std::function<Value const&(typename Samples::value_type const&)> const&
        get_value,
    std::function<Derivative<Value, Argument> const&(
        typename Samples::value_type const&)> const& get_derivative,
    typename Hilbert<Difference<Value>>::NormType const& tolerance);

}  // namespace internal_fit_hermite_spline

using internal_fit_hermite_spline::FitHermiteSpline;

}  // namespace numerics
}  // namespace principia

#include "numerics/fit_hermite_spline_body.hpp"
