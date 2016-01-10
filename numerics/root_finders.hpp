#pragma once

namespace principia {
namespace numerics {

// Approximates a root of |f| between |lower_bound| and |upper_bound| by
// bisection.  The result is less than one ULP from a root of any continuous
// function agreeing with |f| on the values of |Argument|.
// |f(lower_bound)| and |f(upper_bound)| must be nonzero and of opposite signs.
template<typename Argument, typename Function>
Argument Bisect(Function f,
                Argument const& lower_bound,
                Argument const& upper_bound);

}  // namespace numerics
}  // namespace principia

#include "numerics/root_finders_body.hpp"
