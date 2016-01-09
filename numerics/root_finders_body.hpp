#pragma once

#include "root_finders.hpp"

#include "geometry/barycentre_calculator.hpp"
#include "geometry/sign.hpp"
#include "glog/logging.h"

namespace principia {

using geometry::Barycentre;
using geometry::Sign;

namespace numerics {

template<typename Argument, typename Function>
Argument Bisect(Function f,
                Argument const& lower_bound,
                Argument const& upper_bound) {
  using Value = decltype(f(lower_bound));
  Value const zero{};
  Value f_upper = f(upper_bound);
  Value f_lower = f(lower_bound);
  CHECK(f_lower > zero &&  zero > f_upper || f_lower < zero &&  zero < f_upper);
  Argument lower = lower_bound;
  Argument upper = upper_bound;
  for (;;) {
    Argument const middle = Barycentre({lower, upper}, {1, 1});
    // The size of the interval has reached one ULP.
    if (middle == lower || middle == upper) {
      return middle;
    }
    Value const f_middle = f(middle);
    if (f_middle == 0) {
      return middle;
    } else if (Sign(f_middle) == Sign(f_upper)) {
      upper = middle;
    } else {
      lower = middle;
    }
  }
}

}  // namespace numerics
}  // namespace prinicpia
