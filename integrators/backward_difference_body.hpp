
#pragma once

#include "integrators/backward_difference.hpp"

namespace principia {
namespace integrators {

// The coefficients below were computed with Mathematica using the algorithm in
// Fornberg (1988), Generation of Finite Difference Formulas on Arbitrarily
// Spaced Grids.

template<>
inline BackwardDifference<0> const& FirstDerivativeBackwardDifference<0>() {
  static BackwardDifference<0> const backward_difference{{1.0, -1.0}, 1.0};
  return backward_difference;
}

template<>
inline BackwardDifference<1> const& FirstDerivativeBackwardDifference<1>() {
  static BackwardDifference<1> const backward_difference{{3.0, -4.0, 1.0}, 2.0};
  return backward_difference;
}

template<>
inline BackwardDifference<2> const& FirstDerivativeBackwardDifference<2>() {
  static BackwardDifference<2> const backward_difference{
      {11.0, -18.0, 9.0, -2.0}, 6.0};
  return backward_difference;
}

template<>
inline BackwardDifference<3> const& FirstDerivativeBackwardDifference<3>() {
  static BackwardDifference<3> const backward_difference{
      {25.0, -48.0, 36.0, -16.0, 3.0}, 12.0};
  return backward_difference;
}

template<>
inline BackwardDifference<4> const& FirstDerivativeBackwardDifference<4>() {
  static BackwardDifference<4> const backward_difference{
      {137.0, -300.0, 300.0, -200.0, 75.0, -12.0}, 60.0};
  return backward_difference;
}

template<>
inline BackwardDifference<5> const& FirstDerivativeBackwardDifference<5>() {
  static BackwardDifference<5> const backward_difference{
      {147.0, -360.0, 450.0, -400.0, 225.0, -72.0, 10.0}, 60.0};
  return backward_difference;
}

template<>
inline BackwardDifference<6> const& FirstDerivativeBackwardDifference<6>() {
  static BackwardDifference<6> const backward_difference{
      {1089.0, -2940.0, 4410.0, -4900.0, 3675.0, -1764.0, 490.0, -60.0}, 420.0};
  return backward_difference;
}

template<>
inline BackwardDifference<7> const& FirstDerivativeBackwardDifference<7>() {
  static BackwardDifference<7> const backward_difference{{2283.0,
                                                          -6720.0,
                                                          11760.0,
                                                          -15680.0,
                                                          14700.0,
                                                          -9408.0,
                                                          3920.0,
                                                          -960.0,
                                                          105.0},
                                                         840.0};
  return backward_difference;
}

template<>
inline BackwardDifference<8> const& FirstDerivativeBackwardDifference<8>() {
  static BackwardDifference<8> const backward_difference{{7129.0,
                                                          -22680.0,
                                                          45360.0,
                                                          -70560.0,
                                                          79380.0,
                                                          -63504.0,
                                                          35280.0,
                                                          -12960.0,
                                                          2835.0,
                                                          -280.0},
                                                         2520.0};
  return backward_difference;
}

template<>
inline BackwardDifference<9> const& FirstDerivativeBackwardDifference<9>() {
  static BackwardDifference<9> const backward_difference{{7381.0,
                                                          -25200.0,
                                                          56700.0,
                                                          -100800.0,
                                                          132300.0,
                                                          -127008.0,
                                                          88200.0,
                                                          -43200.0,
                                                          14175.0,
                                                          -2800.0,
                                                          252.0},
                                                         2520.0};
  return backward_difference;
}

template<>
inline BackwardDifference<10> const& FirstDerivativeBackwardDifference<10>() {
  static BackwardDifference<10> const backward_difference{{83711.0,
                                                           -304920.0,
                                                           762300.0,
                                                           -1524600.0,
                                                           2286900.0,
                                                           -2561328.0,
                                                           2134440.0,
                                                           -1306800.0,
                                                           571725.0,
                                                           -169400.0,
                                                           30492.0,
                                                           -2520.0},
                                                          27720.0};
  return backward_difference;
}

template<>
inline BackwardDifference<11> const& FirstDerivativeBackwardDifference<11>() {
  static BackwardDifference<11> const backward_difference{{86021.0,
                                                           -332640.0,
                                                           914760.0,
                                                           -2032800.0,
                                                           3430350.0,
                                                           -4390848.0,
                                                           4268880.0,
                                                           -3136320.0,
                                                           1715175.0,
                                                           -677600.0,
                                                           182952.0,
                                                           -30240.0,
                                                           2310.0},
                                                          27720.0};
  return backward_difference;
}

template<>
inline BackwardDifference<12> const& FirstDerivativeBackwardDifference<12>() {
  static BackwardDifference<12> const backward_difference{{1145993.0,
                                                           -4684680.0,
                                                           14054040.0,
                                                           -34354320.0,
                                                           64414350.0,
                                                           -92756664.0,
                                                           103062960.0,
                                                           -88339680.0,
                                                           57972915.0,
                                                           -28628600.0,
                                                           10306296.0,
                                                           -2555280.0,
                                                           390390.0,
                                                           -27720.0},
                                                          360360.0};
  return backward_difference;
}

template<>
inline BackwardDifference<13> const& FirstDerivativeBackwardDifference<13>() {
  static BackwardDifference<13> const backward_difference{{1171733.0,
                                                           -5045040.0,
                                                           16396380.0,
                                                           -43723680.0,
                                                           90180090.0,
                                                           -144288144.0,
                                                           180360180.0,
                                                           -176679360.0,
                                                           135270135.0,
                                                           -80160080.0,
                                                           36072036.0,
                                                           -11924640.0,
                                                           2732730.0,
                                                           -388080.0,
                                                           25740.0},
                                                          360360.0};
  return backward_difference;
}

}  // namespace integrators
}  // namespace principia
