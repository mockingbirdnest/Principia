#pragma once

#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"

namespace principia {
namespace quantities {
namespace _elementary_functions {
namespace internal {

// Equivalent to |std::fma(x, y, z)|.
template<typename Q1, typename Q2>
Product<Q1, Q2> FusedMultiplyAdd(Q1 const& x,
                                 Q2 const& y,
                                 Product<Q1, Q2> const& z);
template<typename Q1, typename Q2>
Product<Q1, Q2> FusedMultiplySubtract(Q1 const& x,
                                      Q2 const& y,
                                      Product<Q1, Q2> const& z);
template<typename Q1, typename Q2>
Product<Q1, Q2> FusedNegatedMultiplyAdd(Q1 const& x,
                                        Q2 const& y,
                                        Product<Q1, Q2> const& z);
template<typename Q1, typename Q2>
Product<Q1, Q2> FusedNegatedMultiplySubtract(Q1 const& x,
                                             Q2 const& y,
                                             Product<Q1, Q2> const& z);

// Equivalent to |std::abs(x)|.
template<typename Q>
Q Abs(Q const& quantity);

// Returns a value between zero and |modulus|.
template<typename Q>
Q Mod(Q const& argument, Q const& modulus);

// Equivalent to |std::sqrt(x)|.
template<typename Q>
SquareRoot<Q> Sqrt(Q const& x);

// Equivalent to |std::cbrt(x)|.
template<typename Q>
CubeRoot<Q> Cbrt(Q const& x);

// Not equivalent to |std::nextafter(x)|; follows IEEE 754:2008 conventions
// instead of C++ ones.  In particular, |NextUp(-0.0) == NextUp(+0.0)|.
template<typename Q>
constexpr Q NextUp(Q const& x);
template<typename Q>
constexpr Q NextDown(Q const& x);

// Equivalent to |std::pow(x, exponent)| unless -3 ≤ x ≤ 3, in which case
// explicit specialization yields multiplications statically.
template<int exponent, typename Q>
constexpr Exponentiation<Q, exponent> Pow(Q const& x);

double Sin(Angle const& α);
double Cos(Angle const& α);
double Tan(Angle const& α);

Angle ArcSin(double x);
Angle ArcCos(double x);
Angle ArcTan(double x);
Angle ArcTan(double y, double x);
template<typename D>
Angle ArcTan(Quantity<D> const& y, Quantity<D> const& x);

// We consider hyperbolic functions as dealing with quotients of arc length to
// curvature radius in the hyperbolic plane, which are angles. This explains the
// use of "arc" for inverse functions.

double Sinh(Angle const& α);
double Cosh(Angle const& α);
double Tanh(Angle const& α);

Angle ArcSinh(double x);
Angle ArcCosh(double x);
Angle ArcTanh(double x);

// Returns the element of {α + 2nπ | n ∈ ℤ} which is closest to
// |previous_angle|.
Angle UnwindFrom(Angle const& previous_angle, Angle const& α);

}  // namespace internal

using internal::Abs;
using internal::ArcCos;
using internal::ArcCosh;
using internal::ArcSin;
using internal::ArcSinh;
using internal::ArcTan;
using internal::ArcTanh;
using internal::Cbrt;
using internal::Cos;
using internal::Cosh;
using internal::FusedMultiplyAdd;
using internal::FusedMultiplySubtract;
using internal::FusedNegatedMultiplyAdd;
using internal::FusedNegatedMultiplySubtract;
using internal::Mod;
using internal::NextDown;
using internal::NextUp;
using internal::Pow;
using internal::Sin;
using internal::Sinh;
using internal::Sqrt;
using internal::Tan;
using internal::Tanh;
using internal::UnwindFrom;

}  // namespace _elementary_functions
}  // namespace quantities
}  // namespace principia

namespace principia::quantities {
using namespace principia::quantities::_elementary_functions;
}  // namespace principia::quantities

#include "quantities/elementary_functions_body.hpp"
