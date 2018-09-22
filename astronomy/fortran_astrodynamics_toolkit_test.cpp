#include "astronomy/fortran_astrodynamics_toolkit.hpp"

namespace principia {
namespace astronomy {
namespace fortran_astrodynamics_toolkit {

auto x = Grav<2, 2>(R3Element<double>(),
                    0.0,
                    0.0,
                    FixedMatrix<double, 2, 3>(),
                    FixedMatrix<double, 2, 3>());

}  // namespace fortran_astrodynamics_toolkit
}  // namespace astronomy
}  // namespace principia
