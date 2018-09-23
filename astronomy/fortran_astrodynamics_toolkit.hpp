
#pragma once

#include "geometry/r3_element.hpp"
#include "numerics/fixed_arrays.hpp"

namespace principia {
namespace astronomy {
namespace fortran_astrodynamics_toolkit {

using geometry::R3Element;
using numerics::FixedMatrix;

template<int nmodel, int mmodel>
R3Element<double> Grav(R3Element<double> const& rgr,
                       double mu,
                       double rbar,
                       FixedMatrix<double, nmodel + 1, nmodel + 1> const& cnm,
                       FixedMatrix<double, nmodel + 1, nmodel + 1> const& snm);

}  // namespace fortran_astrodynamics_toolkit
}  // namespace astronomy
}  // namespace principia

#include "astronomy/fortran_astrodynamics_toolkit_body.hpp"
