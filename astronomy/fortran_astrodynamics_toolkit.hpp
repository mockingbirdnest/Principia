#pragma once

// This code is a straightforward translation of Jacob Williams' implementation
// of [Lea86] and [Lea87].  The original comes with the following notice:
/*
 Fortran Astrodynamics Toolkit
 https://github.com/jacobwilliams/Fortran-Astrodynamics-Toolkit

 Copyright (c) 2014-2018, Jacob Williams
 All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright notice, this
  list of conditions and the following disclaimer in the documentation and/or
  other materials provided with the distribution.

* The names of its contributors may not be used to endorse or promote products
  derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include "geometry/r3_element.hpp"
#include "numerics/fixed_arrays.hpp"

namespace principia {
namespace astronomy {
namespace _fortran_astrodynamics_toolkit {
namespace internal {

using namespace principia::geometry::_r3_element;
using namespace principia::numerics::_fixed_arrays;

template<int nmodel, int mmodel>
R3Element<double> ComputeGravityAccelerationLear(
    R3Element<double> const& rgr,
    double mu,
    double rbar,
    FixedMatrix<double, nmodel + 1, nmodel + 1> const& cnm,
    FixedMatrix<double, nmodel + 1, nmodel + 1> const& snm);

}  // namespace internal

using internal::ComputeGravityAccelerationLear;

}  // namespace _fortran_astrodynamics_toolkit
}  // namespace astronomy
}  // namespace principia

#include "astronomy/fortran_astrodynamics_toolkit_body.hpp"
