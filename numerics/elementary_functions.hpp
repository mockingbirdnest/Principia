#pragma once

#include "base/concepts.hpp"
#include "numerics/fma.hpp"
#include "numerics/m128d.hpp"
#include "quantities/arithmetic.hpp"
#include "quantities/concepts.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"

namespace principia {
namespace numerics {
namespace _elementary_functions {
namespace internal {

using namespace principia::base::_concepts;
using namespace principia::numerics::_fma;
using namespace principia::numerics::_m128d;
using namespace principia::quantities::_arithmetic;
using namespace principia::quantities::_concepts;
using namespace principia::quantities::_named_quantities;
using namespace principia::quantities::_quantities;

// Equivalent to `std::fma(x, y, z)`.
template<FMAPresence fma_presence = FMAPresence::Unknown,
         boost_cpp_int Q1,
         boost_cpp_int Q2>
Product<Q1, Q2> FusedMultiplyAdd(Q1 const& x,
                                 Q2 const& y,
                                 Product<Q1, Q2> const& z);
template<FMAPresence fma_presence = FMAPresence::Unknown,
         boost_cpp_rational Q1,
         boost_cpp_rational Q2>
Product<Q1, Q2> FusedMultiplyAdd(Q1 const& x,
                                 Q2 const& y,
                                 Product<Q1, Q2> const& z);
template<FMAPresence fma_presence = FMAPresence::Unknown,
         boost_cpp_bin_float Q1,
         boost_cpp_bin_float Q2>
Product<Q1, Q2> FusedMultiplyAdd(Q1 const& x,
                                 Q2 const& y,
                                 Product<Q1, Q2> const& z);
template<FMAPresence fma_presence = FMAPresence::Unknown,
         quantity Q1,
         quantity Q2>
Product<Q1, Q2> FusedMultiplyAdd(Q1 const& x,
                                 Q2 const& y,
                                 Product<Q1, Q2> const& z);
template<FMAPresence fma_presence = FMAPresence::Unknown>
M128D FusedMultiplyAdd(M128D const& x, M128D const& y, M128D const& z);

template<FMAPresence fma_presence = FMAPresence::Unknown,
         boost_cpp_int Q1,
         boost_cpp_int Q2>
Product<Q1, Q2> FusedMultiplySubtract(Q1 const& x,
                                      Q2 const& y,
                                      Product<Q1, Q2> const& z);
template<FMAPresence fma_presence = FMAPresence::Unknown,
         boost_cpp_rational Q1,
         boost_cpp_rational Q2>
Product<Q1, Q2> FusedMultiplySubtract(Q1 const& x,
                                      Q2 const& y,
                                      Product<Q1, Q2> const& z);
template<FMAPresence fma_presence = FMAPresence::Unknown,
         boost_cpp_bin_float Q1,
         boost_cpp_bin_float Q2>
Product<Q1, Q2> FusedMultiplySubtract(Q1 const& x,
                                      Q2 const& y,
                                      Product<Q1, Q2> const& z);
template<FMAPresence fma_presence = FMAPresence::Unknown,
         quantity Q1,
         quantity Q2>
Product<Q1, Q2> FusedMultiplySubtract(Q1 const& x,
                                      Q2 const& y,
                                      Product<Q1, Q2> const& z);
template<FMAPresence fma_presence = FMAPresence::Unknown>
M128D FusedMultiplySubtract(M128D const& x, M128D const& y, M128D const& z);

template<FMAPresence fma_presence = FMAPresence::Unknown,
         boost_cpp_int Q1,
         boost_cpp_int Q2>
Product<Q1, Q2> FusedNegatedMultiplyAdd(Q1 const& x,
                                        Q2 const& y,
                                        Product<Q1, Q2> const& z);
template<FMAPresence fma_presence = FMAPresence::Unknown,
         boost_cpp_rational Q1,
         boost_cpp_rational Q2>
Product<Q1, Q2> FusedNegatedMultiplyAdd(Q1 const& x,
                                        Q2 const& y,
                                        Product<Q1, Q2> const& z);
template<FMAPresence fma_presence = FMAPresence::Unknown,
         boost_cpp_bin_float Q1,
         boost_cpp_bin_float Q2>
Product<Q1, Q2> FusedNegatedMultiplyAdd(Q1 const& x,
                                        Q2 const& y,
                                        Product<Q1, Q2> const& z);
template<FMAPresence fma_presence = FMAPresence::Unknown,
         quantity Q1,
         quantity Q2>
Product<Q1, Q2> FusedNegatedMultiplyAdd(Q1 const& x,
                                        Q2 const& y,
                                        Product<Q1, Q2> const& z);
template<FMAPresence fma_presence = FMAPresence::Unknown>
M128D FusedNegatedMultiplyAdd(M128D const& x, M128D const& y, M128D const& z);

template<FMAPresence fma_presence = FMAPresence::Unknown,
         boost_cpp_int Q1,
         boost_cpp_int Q2>
Product<Q1, Q2> FusedNegatedMultiplySubtract(Q1 const& x,
                                             Q2 const& y,
                                             Product<Q1, Q2> const& z);
template<FMAPresence fma_presence = FMAPresence::Unknown,
         boost_cpp_rational Q1,
         boost_cpp_rational Q2>
Product<Q1, Q2> FusedNegatedMultiplySubtract(Q1 const& x,
                                             Q2 const& y,
                                             Product<Q1, Q2> const& z);
template<FMAPresence fma_presence = FMAPresence::Unknown,
         boost_cpp_bin_float Q1,
         boost_cpp_bin_float Q2>
Product<Q1, Q2> FusedNegatedMultiplySubtract(Q1 const& x,
                                             Q2 const& y,
                                             Product<Q1, Q2> const& z);
template<FMAPresence fma_presence = FMAPresence::Unknown,
         quantity Q1,
         quantity Q2>
Product<Q1, Q2> FusedNegatedMultiplySubtract(Q1 const& x,
                                             Q2 const& y,
                                             Product<Q1, Q2> const& z);
template<FMAPresence fma_presence = FMAPresence::Unknown>
M128D FusedNegatedMultiplySubtract(M128D const& x,
                                   M128D const& y,
                                   M128D const& z);


// Equivalent to `std::abs(x)`.
template<typename Q>
Q Abs(Q const& quantity);

// Returns a value between zero and `modulus`.
template<typename Q>
Q Mod(Q const& argument, Q const& modulus);

// Equivalent to `std::sqrt(x)`.
template<typename Q>
SquareRoot<Q> Sqrt(Q const& x);

// Equivalent to `std::cbrt(x)`.
template<typename Q>
CubeRoot<Q> Cbrt(Q const& x);

// Not equivalent to `std::nextafter(x)`; follows IEEE 754:2008 conventions
// instead of C++ ones.  In particular, `NextUp(-0.0) == NextUp(+0.0)`.
template<typename Q>
constexpr Q NextUp(Q const& x);
template<typename Q>
constexpr Q NextDown(Q const& x);

// Equivalent to `std::pow(x, exponent)` unless -3 ≤ x ≤ 3, in which case
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
// `previous_angle`.
Angle UnwindFrom(Angle const& previous_angle, Angle const& α);

// Only dimensionless quantities can be rounded.
template<typename Q>
  requires boost_cpp_number<Q> || std::floating_point<Q>
Q Round(Q const& x);

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
using internal::Round;
using internal::Sin;
using internal::Sinh;
using internal::Sqrt;
using internal::Tan;
using internal::Tanh;
using internal::UnwindFrom;

}  // namespace _elementary_functions
}  // namespace numerics
}  // namespace principia

#include "numerics/elementary_functions_body.hpp"
