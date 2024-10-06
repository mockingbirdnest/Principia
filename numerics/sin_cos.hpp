#pragma once

namespace principia {
namespace numerics {
namespace _sin_cos {
namespace internal {

double Sin(double x);
double Cos(double x);

}  // namespace internal

using internal::Cos;
using internal::Sin;

}  // namespace _sin_cos
}  // namespace numerics
}  // namespace principia
