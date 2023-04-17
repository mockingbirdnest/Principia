#pragma once

namespace principia {
namespace mathematica {
namespace _integrator_plots {
namespace internal {

void GenerateSimpleHarmonicMotionWorkErrorGraphs();
void GenerateKeplerProblemWorkErrorGraphs(double eccentricity);

}  // namespace internal

using internal::GenerateKeplerProblemWorkErrorGraphs;
using internal::GenerateSimpleHarmonicMotionWorkErrorGraphs;

}  // namespace _integrator_plots
}  // namespace mathematica
}  // namespace principia
