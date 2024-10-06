#include "numerics/sin_cos.hpp"

#include <cmath>

#include "numerics/accurate_tables.mathematica.h"
#include "numerics/double_precision.hpp"
#include "numerics/fma.hpp"
#include "numerics/polynomial_evaluators.hpp"

namespace principia {
namespace numerics {
namespace _sin_cos {
namespace internal {

using namespace principia::numerics::_accurate_tables;
using namespace principia::numerics::_double_precision;
using namespace principia::numerics::_fma;
using namespace principia::numerics::_polynomial_evaluators;

using Argument = double;
using Value = double;

constexpr Argument sin_near_zero_cutoff = 1.0 / 1024.0;
constexpr Argument table_spacing_reciprocal = 512.0;

template<FMAPolicy fma_policy>
using Polynomial1 = HornerEvaluator<Value, Argument, 1, fma_policy>;

// TODO(phl): Take the perturbation into account in the polynomials.

template<FMAPolicy fma_policy>
Value SinPolynomial(Argument const x) {
  // Absolute error better than 84.8 bits over an interval of radius 1/1024.
  return Polynomial1<fma_policy>::Evaluate(
      {-0x1.5555555555555'555p-3, 0x1.111110B24ACB5'617p-7}, x);
}

template<FMAPolicy fma_policy>
Value SinPolynomialNearZero(Argument const x) {
  // Relative error better than 74.5 bits over an interval of radius 1/1024.
  return Polynomial1<fma_policy>::Evaluate(
      {-0x1.5555555555555'555p-3, 0x1.111110B40E889'1076p-7}, x);
}

template<FMAPolicy fma_policy>
Value CosPolynomial(Argument const x) {
  // Absolute error better than 72.4 bits over an interval of radius 1/1024.
  return Polynomial1<fma_policy>::Evaluate(
      {-0x1.FFFFFFFFFFFFF'000p-2, 0x1.555554B290E69'14113p-5}, x);
}

template<FMAPolicy fma_policy>
FORCE_INLINE(inline)
Value SinImplementation(Argument const x) {
  if (x < sin_near_zero_cutoff) {
    auto const x² = x * x;
    auto const x³ = x² * x;
    return x + x³ * SinPolynomialNearZero<fma_policy>(x²);
  } else {
    auto const i = std::lround(x * table_spacing_reciprocal);
    auto const& accurate_values = SinCosAccurateTable[i];
    auto const& x₀ = accurate_values.x;
    auto const& sin_x₀ = accurate_values.sin_x;
    auto const& cos_x₀ = accurate_values.cos_x;
    auto const h = x - x₀;

    auto const sin_x₀_plus_h_cos_x₀ =
        TwoProductAdd<fma_policy>(cos_x₀, h, sin_x₀);
    auto const h² = h * h;
    auto const h³ = h² * h;
    return sin_x₀_plus_h_cos_x₀.value +
           ((sin_x₀ * h² * CosPolynomial<fma_policy>(h²) +
             cos_x₀ * h³ * SinPolynomial<fma_policy>(h²)) +
            sin_x₀_plus_h_cos_x₀.error);
  }
}

template<FMAPolicy fma_policy>
FORCE_INLINE(inline)
Value CosImplementation(Argument const x) {
  auto const i = std::lround(x * table_spacing_reciprocal);
  auto const& accurate_values = SinCosAccurateTable[i];
  auto const& x₀ = accurate_values.x;
  auto const& sin_x₀ = accurate_values.sin_x;
  auto const& cos_x₀ = accurate_values.cos_x;
  auto const h = x - x₀;

  auto const cos_x₀_minus_h_sin_x₀ =
      TwoProductNegatedAdd<fma_policy>(sin_x₀, h, cos_x₀);
  auto const h² = h * h;
  auto const h³ = h² * h;
  return cos_x₀_minus_h_sin_x₀.value +
         ((cos_x₀ * h² * CosPolynomial<fma_policy>(h²) -
           sin_x₀ * h³ * SinPolynomial<fma_policy>(h²)) +
          cos_x₀_minus_h_sin_x₀.error);
}

Value Sin(Argument const x) {
  return UseHardwareFMA ? SinImplementation<FMAPolicy::Force>(x)
                        : SinImplementation<FMAPolicy::Disallow>(x);
}

Value Cos(Argument const x) {
  return UseHardwareFMA ? CosImplementation<FMAPolicy::Force>(x)
                        : CosImplementation<FMAPolicy::Disallow>(x);
}

}  // namespace internal
}  // namespace _sin_cos
}  // namespace numerics
}  // namespace principia
