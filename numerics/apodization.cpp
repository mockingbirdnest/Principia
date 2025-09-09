#include "numerics/apodization.hpp"

#include "geometry/barycentre_calculator.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/si.hpp"

// The formulæ below differs from those in
// https://en.wikipedia.org/wiki/Window_function because our functions have as
// their origin the centre of the interval, while those in Wikipedia have as
// their origin the lower bound of the interval.  Specifically consider:
//   f(t) = a₀ - a₁ cos ωt - a₂ cos 2ωt - a₃ cos 3ωt - a₄ cos 4ωt
// for t ∈ [0, 2π/ω].  Rewriting this formula in terms of τ = t - π/ω yields:
//   f(τ) = a₀ + a₁ cos ωτ - a₂ cos 2ωτ + a₃ cos 3ωτ - a₄ cos 4ωτ
// for τ ∈ [-π/ω, π/ω].

namespace principia {
namespace numerics {
namespace _apodization {
namespace internal {

using namespace principia::geometry::_barycentre_calculator;
using namespace principia::quantities::_named_quantities;
using namespace principia::quantities::_si;

PoissonSeries<double, 0, 0> Dirichlet(Instant const& t_min,
                                      Instant const& t_max) {
  using Result = PoissonSeries<double, 0, 0>;
  using AperiodicPolynomial = typename Result::AperiodicPolynomial;
  Instant const t_mid = Barycentre({t_min, t_max});
  return Result(AperiodicPolynomial({1}, t_mid), {});
}

PoissonSeries<double, 0, 0> Sine(Instant const& t_min,
                                 Instant const& t_max) {
  using Result = PoissonSeries<double, 0, 0>;
  using AperiodicPolynomial = typename Result::AperiodicPolynomial;
  using PeriodicPolynomial = typename Result::PeriodicPolynomial;
  AngularFrequency const ω = π * Radian / (t_max - t_min);
  Instant const t_mid = Barycentre({t_min, t_max});
  return Result(AperiodicPolynomial({0}, t_mid),
                {{ω,
                  {.sin = PeriodicPolynomial({0}, t_mid),
                   .cos = PeriodicPolynomial({1}, t_mid)}}});
}

PoissonSeries<double, 0, 0> Hann(Instant const& t_min,
                                 Instant const& t_max) {
  using Result = PoissonSeries<double, 0, 0>;
  using AperiodicPolynomial = typename Result::AperiodicPolynomial;
  using PeriodicPolynomial = typename Result::PeriodicPolynomial;
  AngularFrequency const ω = 2 * π * Radian / (t_max - t_min);
  Instant const t_mid = Barycentre({t_min, t_max});
  return Result(AperiodicPolynomial({0.5}, t_mid),
                {{ω,
                  {.sin = PeriodicPolynomial({0}, t_mid),
                   .cos = PeriodicPolynomial({0.5}, t_mid)}}});
}

PoissonSeries<double, 0, 0> Hamming(Instant const& t_min,
                                    Instant const& t_max) {
  using Result = PoissonSeries<double, 0, 0>;
  AngularFrequency const ω = 2 * π * Radian / (t_max - t_min);
  using AperiodicPolynomial = typename Result::AperiodicPolynomial;
  using PeriodicPolynomial = typename Result::PeriodicPolynomial;
  Instant const t_mid = Barycentre({t_min, t_max});
  return Result(AperiodicPolynomial({25.0 / 46.0}, t_mid),
                {{ω,
                  {.sin = PeriodicPolynomial({0}, t_mid),
                   .cos = PeriodicPolynomial({21.0 / 46.0}, t_mid)}}});
}

PoissonSeries<double, 0, 0> Blackman(Instant const& t_min,
                                     Instant const& t_max) {
  using Result = PoissonSeries<double, 0, 0>;
  using AperiodicPolynomial = typename Result::AperiodicPolynomial;
  using PeriodicPolynomial = typename Result::PeriodicPolynomial;
  AngularFrequency const ω = 2 * π * Radian / (t_max - t_min);
  Instant const t_mid = Barycentre({t_min, t_max});
  return Result(AperiodicPolynomial({0.42}, t_mid),
                {{ω,
                  {.sin = PeriodicPolynomial({0}, t_mid),
                   .cos = PeriodicPolynomial({0.5}, t_mid)}},
                 {2 * ω,
                  {.sin = PeriodicPolynomial({0}, t_mid),
                   .cos = PeriodicPolynomial({0.08}, t_mid)}}});
}

PoissonSeries<double, 0, 0> ExactBlackman(Instant const& t_min,
                                          Instant const& t_max) {
  using Result = PoissonSeries<double, 0, 0>;
  using AperiodicPolynomial = typename Result::AperiodicPolynomial;
  using PeriodicPolynomial = typename Result::PeriodicPolynomial;
  AngularFrequency const ω = 2 * π * Radian / (t_max - t_min);
  Instant const t_mid = Barycentre({t_min, t_max});
  return Result(AperiodicPolynomial({3969.0 / 9304.0}, t_mid),
                {{ω,
                  {.sin = PeriodicPolynomial({0}, t_mid),
                   .cos = PeriodicPolynomial({1155.0 / 2326.0}, t_mid)}},
                 {2 * ω,
                  {.sin = PeriodicPolynomial({0}, t_mid),
                   .cos = PeriodicPolynomial({715.0 / 9304.0}, t_mid)}}});
}

PoissonSeries<double, 0, 0> Nuttall(Instant const& t_min,
                                    Instant const& t_max) {
  using Result = PoissonSeries<double, 0, 0>;
  using AperiodicPolynomial = typename Result::AperiodicPolynomial;
  using PeriodicPolynomial = typename Result::PeriodicPolynomial;
  AngularFrequency const ω = 2 * π * Radian / (t_max - t_min);
  Instant const t_mid = Barycentre({t_min, t_max});
  return Result(AperiodicPolynomial({0.355768}, t_mid),
                {{ω,
                  {.sin = PeriodicPolynomial({0}, t_mid),
                   .cos = PeriodicPolynomial({0.487396}, t_mid)}},
                 {2 * ω,
                  {.sin = PeriodicPolynomial({0}, t_mid),
                   .cos = PeriodicPolynomial({0.144232}, t_mid)}},
                 {3 * ω,
                  {.sin = PeriodicPolynomial({0}, t_mid),
                   .cos = PeriodicPolynomial({0.012604}, t_mid)}}});
}

PoissonSeries<double, 0, 0> BlackmanNuttall(Instant const& t_min,
                                            Instant const& t_max) {
  using Result = PoissonSeries<double, 0, 0>;
  using AperiodicPolynomial = typename Result::AperiodicPolynomial;
  using PeriodicPolynomial = typename Result::PeriodicPolynomial;
  AngularFrequency const ω = 2 * π * Radian / (t_max - t_min);
  Instant const t_mid = Barycentre({t_min, t_max});
  return Result(AperiodicPolynomial({0.3635819}, t_mid),
                {{ω,
                  {.sin = PeriodicPolynomial({0}, t_mid),
                   .cos = PeriodicPolynomial({0.4891775}, t_mid)}},
                 {2 * ω,
                  {.sin = PeriodicPolynomial({0}, t_mid),
                   .cos = PeriodicPolynomial({0.1365995}, t_mid)}},
                 {3 * ω,
                  {.sin = PeriodicPolynomial({0}, t_mid),
                   .cos = PeriodicPolynomial({0.0106411}, t_mid)}}});
}

PoissonSeries<double, 0, 0> BlackmanHarris(Instant const& t_min,
                                           Instant const& t_max) {
  using Result = PoissonSeries<double, 0, 0>;
  using AperiodicPolynomial = typename Result::AperiodicPolynomial;
  using PeriodicPolynomial = typename Result::PeriodicPolynomial;
  AngularFrequency const ω = 2 * π * Radian / (t_max - t_min);
  Instant const t_mid = Barycentre({t_min, t_max});
  return Result(AperiodicPolynomial({0.35875}, t_mid),
                {{ω,
                  {.sin = PeriodicPolynomial({0}, t_mid),
                   .cos = PeriodicPolynomial({0.48829}, t_mid)}},
                 {2 * ω,
                  {.sin = PeriodicPolynomial({0}, t_mid),
                   .cos = PeriodicPolynomial({0.14128}, t_mid)}},
                 {3 * ω,
                  {.sin = PeriodicPolynomial({0}, t_mid),
                   .cos = PeriodicPolynomial({0.01168}, t_mid)}}});
}

PoissonSeries<double, 0, 0> ISO18431_2(Instant const& t_min,
                                       Instant const& t_max) {
  using Result = PoissonSeries<double, 0, 0>;
  using AperiodicPolynomial = typename Result::AperiodicPolynomial;
  using PeriodicPolynomial = typename Result::PeriodicPolynomial;
  AngularFrequency const ω = 2 * π * Radian / (t_max - t_min);
  Instant const t_mid = Barycentre({t_min, t_max});
  return Result(AperiodicPolynomial({1.0 / 4.63867187}, t_mid),
                {{ω,
                  {.sin = PeriodicPolynomial({0}, t_mid),
                   .cos = PeriodicPolynomial(
                       {1.93261719 / 4.63867187}, t_mid)}},
                 {2 * ω,
                  {.sin = PeriodicPolynomial({0}, t_mid),
                   .cos = PeriodicPolynomial(
                       {1.28613281 / 4.63867187}, t_mid)}},
                 {3 * ω,
                  {.sin = PeriodicPolynomial({0}, t_mid),
                   .cos = PeriodicPolynomial(
                       {0.38769531 / 4.63867187}, t_mid)}},
                 {4 * ω,
                  {.sin = PeriodicPolynomial({0}, t_mid),
                   .cos = PeriodicPolynomial(
                       {0.03222656 / 4.63867187}, t_mid)}}});
}

}  // namespace internal
}  // namespace _apodization
}  // namespace numerics
}  // namespace principia
