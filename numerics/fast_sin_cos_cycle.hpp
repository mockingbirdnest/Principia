
#pragma once

namespace principia {
namespace numerics {

// An implementation of sin and cos that favours performance at the expense of
// accuracy.  It has an absolute error of roughly 6e-7 and takes about 25
// cycles.  Note that the argument x is not an angle but a number of cycles.
void FastSinCosCycle(double x, double& sin, double& cos);

}  // namespace numerics
}  // namespace principia
