#pragma once

#include "numerics/poisson_series.hpp"

#include <algorithm>
#include <functional>
#include <map>
#include <optional>
#include <utility>
#include <vector>

#include "numerics/double_precision.hpp"
#include "numerics/elementary_functions.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace numerics {
namespace _poisson_series {
namespace internal {

using namespace principia::numerics::_double_precision;
using namespace principia::numerics::_elementary_functions;
using namespace principia::quantities::_si;

// These parameters have been tuned for approximation of the Moon over 3 months
// with 10 periods.

// In the fast/slow algorithms, specifies the maximum number of periods over the
// time interval below which we use Clenshaw-Curtis integration.
constexpr int clenshaw_curtis_max_periods_overall = 40;

// The minimum value of the max_point parameter passed to Clenshaw-Curtis
// integration, irrespective of the frequencies of the argument function.
constexpr int clenshaw_curtis_min_points_overall = 65;

// The maximum number of points use in Clenshaw-Curtis integration for each
// period of the highest frequency of the argument function.
constexpr int clenshaw_curtis_points_per_period = 4;

// The desired relative error on Clenshaw-Curtis integration, as determined by
// two successive computations with increasing number of points.
constexpr double clenshaw_curtis_relative_error = 0x1p-32;

// This function computes ∫ₜ₁ᵗ²(p(t) sin ω t + q(t) cos ω t) dt.
template<typename Value, int degree>
Primitive<Value, Time> AngularFrequencyIntegrate(
    AngularFrequency const& ω,
    PolynomialInMonomialBasis<Value, Instant, degree> const& p,
    PolynomialInMonomialBasis<Value, Instant, degree> const& q,
    Instant const& t1,
    Instant const& t2,
    double const sin_ωt1,
    double const cos_ωt1,
    double const sin_ωt2,
    double const cos_ωt2) {
  static_assert(degree >= 0);
  DoublePrecision<Value> sum;
  sum += q(t2) * sin_ωt2;
  sum -= p(t2) * cos_ωt2;
  sum -= q(t1) * sin_ωt1;
  sum += p(t1) * cos_ωt1;
  if constexpr (degree > 0) {
    sum += AngularFrequencyIntegrate(ω,
                                     /*p=*/-q.template Derivative<1>(),
                                     /*q=*/p.template Derivative<1>(),
                                     t1, t2,
                                     sin_ωt1, cos_ωt1,
                                     sin_ωt2, cos_ωt2);
  }
  return (sum.value + sum.error) / ω * Radian;
}

// This function computes ∫(p(t) sin ω t + q(t) cos ω t) dt where p and q are
// the two parts of the polynomials argument.
template<typename Value,
         int aperiodic_degree, int periodic_degree>
typename PoissonSeries<Primitive<Value, Time>,
                       aperiodic_degree + 1, periodic_degree + 1>::Polynomials
AngularFrequencyPrimitive(
    AngularFrequency const& ω,
    typename PoissonSeries<Value, aperiodic_degree, periodic_degree>::
        Polynomials const& polynomials) {
  using Result = PoissonSeries<Primitive<Value, Time>,
                               aperiodic_degree + 1, periodic_degree + 1>;
  using PeriodicPolynomial = typename Result::PeriodicPolynomial;

  // Integration by parts.
  typename Result::Polynomials const first_part{
      .sin = PeriodicPolynomial(polynomials.cos / ω * Radian),
      .cos = PeriodicPolynomial(-polynomials.sin / ω * Radian)};
  if constexpr (periodic_degree == 0) {
    return first_part;
  } else {
    auto const sin_polynomial =
        -polynomials.cos.template Derivative<1>() / ω * Radian;
    auto const cos_polynomial =
        polynomials.sin.template Derivative<1>() / ω * Radian;
    auto const second_part =
        AngularFrequencyPrimitive<Value,
                                  aperiodic_degree - 1, periodic_degree - 1>(
            ω,
            {.sin = sin_polynomial,
            .cos = cos_polynomial});
    return {.sin = first_part.sin + second_part.sin,
            .cos = first_part.cos + second_part.cos};
  }
}

// A helper for multiplication of Poisson series and pointwise inner product.
// The functor Product must take a pair of Poisson series with the types of left
// and right and return a suitable Poisson series.
template<typename LValue, typename RValue,
         int aperiodic_ldegree, int periodic_ldegree,
         int aperiodic_rdegree, int periodic_rdegree,
         typename Product>
auto Multiply(PoissonSeries<LValue,
                            aperiodic_ldegree, periodic_ldegree> const& left,
              PoissonSeries<RValue,
                            aperiodic_rdegree, periodic_rdegree> const& right,
              Product const& product) {
  using Result = PoissonSeries<
      typename std::invoke_result_t<
          Product,
          typename PoissonSeries<LValue,
                                 aperiodic_ldegree, periodic_ldegree>::
              AperiodicPolynomial,
          typename PoissonSeries<RValue,
                                 aperiodic_rdegree, periodic_rdegree>::
              AperiodicPolynomial>::
          Value,
      std::max({aperiodic_ldegree + aperiodic_rdegree,
                aperiodic_ldegree + periodic_rdegree,
                periodic_ldegree + aperiodic_rdegree,
                periodic_ldegree + periodic_rdegree}),
      std::max({aperiodic_ldegree + periodic_rdegree,
                periodic_ldegree + aperiodic_rdegree,
                periodic_ldegree + periodic_rdegree})>;

  auto aperiodic = typename Result::AperiodicPolynomial(
      product(left.aperiodic_, right.aperiodic_));

  // Compute all the individual terms using elementary trigonometric identities
  // and put them in a vector, because the same frequency may appear multiple
  // times.
  typename Result::PolynomialsByAngularFrequency periodic;
  periodic.reserve(left.periodic_.size() + right.periodic_.size() +
                   2 * left.periodic_.size() * right.periodic_.size());
  for (auto const& [ω, polynomials] : left.periodic_) {
    auto const cos_polynomial = typename Result::PeriodicPolynomial(
        product(polynomials.cos, right.aperiodic_));
    auto const sin_polynomial = typename Result::PeriodicPolynomial(
        product(polynomials.sin, right.aperiodic_));
    periodic.emplace_back(ω,
                          typename Result::Polynomials{.sin = sin_polynomial,
                                                       .cos = cos_polynomial});
  }
  for (auto const& [ω, polynomials] : right.periodic_) {
    auto const cos_polynomial = typename Result::PeriodicPolynomial(
        product(left.aperiodic_, polynomials.cos));
    auto const sin_polynomial = typename Result::PeriodicPolynomial(
        product(left.aperiodic_, polynomials.sin));
    periodic.emplace_back(ω,
                          typename Result::Polynomials{.sin = sin_polynomial,
                                                       .cos = cos_polynomial});
  }
  for (auto const& [ωl, polynomials_left] : left.periodic_) {
    for (auto const& [ωr, polynomials_right] : right.periodic_) {
      auto const cos_cos = typename Result::PeriodicPolynomial(
          product(polynomials_left.cos, polynomials_right.cos));
      auto const cos_sin = typename Result::PeriodicPolynomial(
          product(polynomials_left.cos, polynomials_right.sin));
      auto const sin_cos = typename Result::PeriodicPolynomial(
          product(polynomials_left.sin, polynomials_right.cos));
      auto const sin_sin = typename Result::PeriodicPolynomial(
          product(polynomials_left.sin, polynomials_right.sin));
      periodic.emplace_back(
          ωl - ωr,
          typename Result::Polynomials{.sin = (-cos_sin + sin_cos) / 2,
                                       .cos = (sin_sin + cos_cos) / 2});
      periodic.emplace_back(
          ωl + ωr,
          typename Result::Polynomials{.sin = (cos_sin + sin_cos) / 2,
                                       .cos = (-sin_sin + cos_cos) / 2});
    }
  }

  return Result(typename Result::PrivateConstructor{},
                std::move(aperiodic),
                std::move(periodic));
}

template<typename Value,
         int aperiodic_degree_, int periodic_degree_>
PoissonSeries<Value, aperiodic_degree_, periodic_degree_>::
PoissonSeries(AperiodicPolynomial const& aperiodic,
              PolynomialsByAngularFrequency const& periodic)
    : PoissonSeries(PrivateConstructor{},
                    AperiodicPolynomial(aperiodic),
                    PolynomialsByAngularFrequency(periodic)) {}

template<typename Value,
         int aperiodic_degree_, int periodic_degree_>
template<int higher_aperiodic_degree, int higher_periodic_degree>
PoissonSeries<Value, aperiodic_degree_, periodic_degree_>::
operator PoissonSeries<Value, higher_aperiodic_degree, higher_periodic_degree>()
    const {
  static_assert(aperiodic_degree_ <= higher_aperiodic_degree);
  static_assert(periodic_degree_ <= higher_periodic_degree);
  using Result = PoissonSeries<Value,
                               higher_aperiodic_degree, higher_periodic_degree>;
  auto aperiodic = typename Result::AperiodicPolynomial(aperiodic_);
  typename Result::PolynomialsByAngularFrequency periodic;
  periodic.reserve(periodic_.size());
  for (auto const& [ω, polynomials] : periodic_) {
    periodic.emplace_back(
        ω,
        typename Result::Polynomials{
            .sin = typename Result::PeriodicPolynomial(polynomials.sin),
            .cos = typename Result::PeriodicPolynomial(polynomials.cos)});
  }
  return Result(typename Result::TrustedPrivateConstructor{},
                std::move(aperiodic),
                std::move(periodic));
}

template<typename Value,
         int aperiodic_degree_, int periodic_degree_>
template<int aperiodic_rdegree, int periodic_rdegree>
PoissonSeries<Value, aperiodic_degree_, periodic_degree_>&
PoissonSeries<Value, aperiodic_degree_, periodic_degree_>::
operator+=(PoissonSeries<Value,
                         aperiodic_rdegree, periodic_rdegree> const& right) {
  static_assert(aperiodic_rdegree <= aperiodic_degree_);
  static_assert(periodic_rdegree <= periodic_degree_);
  *this = *this + right;
  return *this;
}

template<typename Value,
         int aperiodic_degree_, int periodic_degree_>
template<int aperiodic_rdegree, int periodic_rdegree>
PoissonSeries<Value, aperiodic_degree_, periodic_degree_>&
PoissonSeries<Value, aperiodic_degree_, periodic_degree_>::
operator-=(PoissonSeries<Value,
                         aperiodic_rdegree, periodic_rdegree> const& right) {
  static_assert(aperiodic_rdegree <= aperiodic_degree_);
  static_assert(periodic_rdegree <= periodic_degree_);
  *this = *this - right;
  return *this;
}

template<typename Value,
         int aperiodic_degree_, int periodic_degree_>
Instant const&
PoissonSeries<Value, aperiodic_degree_, periodic_degree_>::
origin() const {
  return origin_;
}

template<typename Value,
         int aperiodic_degree_, int periodic_degree_>
AngularFrequency
PoissonSeries<Value, aperiodic_degree_, periodic_degree_>::
max_ω() const {
  return periodic_.empty() ? AngularFrequency{} : periodic_.back().first;
}

template<typename Value,
         int aperiodic_degree_, int periodic_degree_>
Value PoissonSeries<Value, aperiodic_degree_, periodic_degree_>::
operator()(Instant const& t) const {
  Time const Δt = t - origin_;
  Value result = aperiodic_(t);
  for (auto const& [ω, polynomials] : periodic_) {
    Angle const ωΔt = ω * Δt;
    result += polynomials.sin(t) * Sin(ωΔt) + polynomials.cos(t) * Cos(ωΔt);
  }
  return result;
}

template<typename Value,
         int aperiodic_degree_, int periodic_degree_>
PoissonSeries<Value, aperiodic_degree_, periodic_degree_>
PoissonSeries<Value, aperiodic_degree_, periodic_degree_>::AtOrigin(
    Instant const& origin) const {
  Time const shift = origin - origin_;
  auto aperiodic = aperiodic_.AtOrigin(origin);

  PolynomialsByAngularFrequency periodic;
  periodic.reserve(periodic_.size());
  for (auto const& [ω, polynomials] : periodic_) {
    double const cos_ω_shift = Cos(ω * shift);
    double const sin_ω_shift = Sin(ω * shift);
    PeriodicPolynomial const sin_at_origin = polynomials.sin.AtOrigin(origin);
    PeriodicPolynomial const cos_at_origin = polynomials.cos.AtOrigin(origin);
    periodic.emplace_back(ω,
                          Polynomials{.sin = sin_at_origin * cos_ω_shift -
                                             cos_at_origin * sin_ω_shift,
                                      .cos = sin_at_origin * sin_ω_shift +
                                             cos_at_origin * cos_ω_shift});
  }
  return {TrustedPrivateConstructor{},
          std::move(aperiodic),
          std::move(periodic)};
}

template<typename Value,
         int aperiodic_degree_, int periodic_degree_>
PoissonSeries<Primitive<Value, Time>,
              aperiodic_degree_ + 1, periodic_degree_ + 1>
PoissonSeries<Value, aperiodic_degree_, periodic_degree_>::
Primitive() const {
  using Result =
      PoissonSeries<quantities::_named_quantities::Primitive<Value, Time>,
                    aperiodic_degree_ + 1, periodic_degree_ + 1>;
  typename Result::AperiodicPolynomial aperiodic = aperiodic_.Primitive();
  typename Result::PolynomialsByAngularFrequency periodic;
  periodic.reserve(periodic_.size());
  for (auto const& [ω, polynomials] : periodic_) {
    periodic.emplace_back(
        ω,
        AngularFrequencyPrimitive<Value,
                                  aperiodic_degree_, periodic_degree_>(
            ω, polynomials));
  }
  return Result{aperiodic, periodic};
}

template<typename Value,
         int aperiodic_degree_, int periodic_degree_>
Primitive<Value, Time>
PoissonSeries<Value, aperiodic_degree_, periodic_degree_>::
Integrate(Instant const& t1,
          Instant const& t2) const {
  quantities::_named_quantities::Primitive<Value, Time> result =
      aperiodic_.Integrate(t1, t2);
  for (auto const& [ω, polynomials] : periodic_) {
    // This implementation follows [HO09], Theorem 1 and [INO06] equation 4.
    // The trigonometric functions are computed only once as we iterate through
    // the degree of the polynomials.
    auto const sin_ωt1 = Sin(ω * (t1 - origin_));
    auto const cos_ωt1 = Cos(ω * (t1 - origin_));
    auto const sin_ωt2 = Sin(ω * (t2 - origin_));
    auto const cos_ωt2 = Cos(ω * (t2 - origin_));
    result += AngularFrequencyIntegrate(ω,
                                        polynomials.sin, polynomials.cos,
                                        t1, t2,
                                        sin_ωt1, cos_ωt1,
                                        sin_ωt2, cos_ωt2);
  }
  return result;
}

template<typename Value,
         int aperiodic_degree_, int periodic_degree_>
template<int aperiodic_wdegree, int periodic_wdegree>
typename Hilbert<Value>::NormType
PoissonSeries<Value, aperiodic_degree_, periodic_degree_>::
Norm(PoissonSeries<double,
                   aperiodic_wdegree, periodic_wdegree> const& weight,
     Instant const& t_min,
     Instant const& t_max) const {
  AngularFrequency const ω_cutoff =
      2 * π * Radian * clenshaw_curtis_max_periods_overall / (t_max - t_min);
  auto const split = Split(ω_cutoff);

  AngularFrequency const max_ω = 2 * split.slow.max_ω() + weight.max_ω();
  std::optional<int> const max_points =
      _quadrature::MaxPointsHeuristicsForAutomaticClenshawCurtis(
          max_ω,
          t_max - t_min,
          clenshaw_curtis_min_points_overall,
          clenshaw_curtis_points_per_period);

  auto slow_integrand = [&split, &weight](Instant const& t) {
    return Hilbert<Value>::Norm²(split.slow(t)) * weight(t);
  };
  auto const slow_quadrature = _quadrature::AutomaticClenshawCurtis(
      slow_integrand,
      t_min, t_max,
      /*max_relative_error=*/clenshaw_curtis_relative_error,
      /*max_points=*/max_points);

  auto const fast_integrand =
      PointwiseInnerProduct(split.fast, split.fast + 2 * split.slow) *
      weight;
  auto const fast_quadrature = fast_integrand.Integrate(t_min, t_max);

  return Sqrt((slow_quadrature + fast_quadrature) / (t_max - t_min));
}

template<typename Value,
         int aperiodic_degree_, int periodic_degree_>
void PoissonSeries<Value, aperiodic_degree_, periodic_degree_>::
WriteToMessage(not_null<serialization::PoissonSeries*> const message) const {
  aperiodic_.WriteToMessage(message->mutable_aperiodic());
  for (auto const& [ω, polynomials] : periodic_) {
    auto* const polynomials_and_angular_frequency = message->add_periodic();
    ω.WriteToMessage(
        polynomials_and_angular_frequency->mutable_angular_frequency());
    polynomials.sin.WriteToMessage(
        polynomials_and_angular_frequency->mutable_sin());
    polynomials.cos.WriteToMessage(
        polynomials_and_angular_frequency->mutable_cos());
  }
}

template<typename Value,
         int aperiodic_degree_, int periodic_degree_>
PoissonSeries<Value, aperiodic_degree_, periodic_degree_>
PoissonSeries<Value, aperiodic_degree_, periodic_degree_>::
ReadFromMessage(serialization::PoissonSeries const& message) {
  auto const aperiodic =
      AperiodicPolynomial::ReadFromMessage(message.aperiodic());
  PolynomialsByAngularFrequency periodic;
  for (auto const& polynomial_and_angular_frequency : message.periodic()) {
    auto const ω = AngularFrequency::ReadFromMessage(
        polynomial_and_angular_frequency.angular_frequency());
    auto const sin = PeriodicPolynomial::ReadFromMessage(
        polynomial_and_angular_frequency.sin());
    auto const cos = PeriodicPolynomial::ReadFromMessage(
        polynomial_and_angular_frequency.cos());
    periodic.emplace_back(ω, Polynomials{.sin = sin, .cos = cos});
  }
  return PoissonSeries(aperiodic, periodic);
}

template<typename Value,
         int aperiodic_degree_, int periodic_degree_>
PoissonSeries<Value, aperiodic_degree_, periodic_degree_>::
PoissonSeries(PrivateConstructor,
              AperiodicPolynomial aperiodic,
              PolynomialsByAngularFrequency periodic)
    : origin_(aperiodic.origin()),
      aperiodic_(std::move(aperiodic)),
      periodic_(std::move(periodic)) {
  // The `periodic` vector may have elements with positive or negative angular
  // frequencies.  Normalize our member variable to only have positive angular
  // frequencies.

  // Sort by ascending frequency, irrespective of sign.
  std::sort(
      periodic_.begin(),
      periodic_.end(),
      [](typename PolynomialsByAngularFrequency::value_type const& left,
         typename PolynomialsByAngularFrequency::value_type const& right) {
        return Abs(left.first) < Abs(right.first);
      });

  // Group the terms together by frequency, removing consecutive terms with the
  // same frequency, normalizing negative frequencies, and moving zero
  // frequencies to the aperiodic term.
  std::optional<AngularFrequency> previous_abs_ω;
  for (auto it = periodic_.begin(); it != periodic_.end();) {
    auto& ω = it->first;
    auto const abs_ω = Abs(ω);
    auto& polynomials = it->second;

    // All polynomials must have the same origin.
    CHECK_EQ(origin_, polynomials.sin.origin());
    CHECK_EQ(origin_, polynomials.cos.origin());

    if (ω < AngularFrequency{}) {
      if (previous_abs_ω.has_value() && previous_abs_ω.value() == -ω) {
        auto& previous_polynomials = std::prev(it)->second;
        previous_polynomials.sin -= polynomials.sin;
        previous_polynomials.cos += polynomials.cos;
        it = periodic_.erase(it);
      } else {
        ω = -ω;
        polynomials.sin = -polynomials.sin;
        ++it;
      }
    } else if (ω > AngularFrequency{}) {
      if (previous_abs_ω.has_value() && previous_abs_ω.value() == ω) {
        auto& previous_polynomials = std::prev(it)->second;
        previous_polynomials.sin += polynomials.sin;
        previous_polynomials.cos += polynomials.cos;
        it = periodic_.erase(it);
      } else {
        ++it;
      }
    } else if constexpr (aperiodic_degree_ >= periodic_degree_) {
      aperiodic_ += AperiodicPolynomial(polynomials.cos);
      it = periodic_.erase(it);
    } else {
      LOG(FATAL) << "Degrees mismatch for zero frequency: " << polynomials.cos;
    }
    previous_abs_ω = abs_ω;
  }
}

template<typename Value,
         int aperiodic_degree_, int periodic_degree_>
PoissonSeries<Value, aperiodic_degree_, periodic_degree_>::
PoissonSeries(TrustedPrivateConstructor,
              AperiodicPolynomial aperiodic,
              PolynomialsByAngularFrequency periodic)
    : origin_(aperiodic.origin()),
      aperiodic_(std::move(aperiodic)),
      periodic_(std::move(periodic)) {}

template<typename Value,
         int aperiodic_degree_, int periodic_degree_>
typename PoissonSeries<Value,
                       aperiodic_degree_, periodic_degree_>::SplitPoissonSeries
PoissonSeries<Value, aperiodic_degree_, periodic_degree_>::
Split(AngularFrequency const& ω_cutoff) const {
  // TODO(phl): Should we try to avoid a linear search and copies?
  typename PoissonSeries::PolynomialsByAngularFrequency slow_periodic;
  typename PoissonSeries::PolynomialsByAngularFrequency fast_periodic;
  for (auto const& [ω, polynomials] : periodic_) {
    if (ω <= ω_cutoff) {
      slow_periodic.emplace_back(ω, polynomials);
    } else {
      fast_periodic.emplace_back(ω, polynomials);
    }
  }

  // The reason for having the slow/fast split is to handle cancellations that
  // occurs between the aperiodic component and the components with low
  // frequencies.  If there are no low frequencies, we might as well put the
  // aperiodic component in the high-frequencies half.
  if (slow_periodic.empty()) {
    PoissonSeries slow(TrustedPrivateConstructor{},
                       typename PoissonSeries::AperiodicPolynomial(
                           {}, aperiodic_.origin()),
                       std::move(slow_periodic));
    PoissonSeries fast(TrustedPrivateConstructor{},
                       aperiodic_,
                       std::move(fast_periodic));
    return {.slow = std::move(slow), .fast = std::move(fast)};
  } else {
    PoissonSeries slow(TrustedPrivateConstructor{},
                       aperiodic_,
                       std::move(slow_periodic));
    PoissonSeries fast(TrustedPrivateConstructor{},
                       typename PoissonSeries::AperiodicPolynomial(
                           {}, aperiodic_.origin()),
                       std::move(fast_periodic));
    return {.slow = std::move(slow), .fast = std::move(fast)};
  }
}

template<typename Value,
         int aperiodic_rdegree, int periodic_rdegree>
PoissonSeries<Value, aperiodic_rdegree, periodic_rdegree>
operator+(PoissonSeries<Value,
                        aperiodic_rdegree, periodic_rdegree> const& right) {
  return right;
}

template<typename Value,
         int aperiodic_rdegree, int periodic_rdegree>
PoissonSeries<Value, aperiodic_rdegree, periodic_rdegree>
operator-(PoissonSeries<Value,
                        aperiodic_rdegree, periodic_rdegree> const& right) {
  using Result = PoissonSeries<Value,
                               aperiodic_rdegree, periodic_rdegree>;
  auto aperiodic = -right.aperiodic_;
  typename Result::PolynomialsByAngularFrequency periodic;
  periodic.reserve(right.periodic_.size());
  for (auto const& [ω, polynomials] : right.periodic_) {
    periodic.emplace_back(
        ω,
        typename Result::Polynomials{.sin = -polynomials.sin,
                                     .cos = -polynomials.cos});
  }
  return {typename Result::TrustedPrivateConstructor{},
          std::move(aperiodic),
          std::move(periodic)};
}

template<typename Value,
         int aperiodic_ldegree, int periodic_ldegree,
         int aperiodic_rdegree, int periodic_rdegree>
PoissonSeries<Value,
              std::max(aperiodic_ldegree, aperiodic_rdegree),
              std::max(periodic_ldegree, periodic_rdegree)>
operator+(PoissonSeries<Value,
                        aperiodic_ldegree, periodic_ldegree> const& left,
          PoissonSeries<Value,
                        aperiodic_rdegree, periodic_rdegree> const& right) {
  using Result =
      PoissonSeries<Value,
                    std::max(aperiodic_ldegree, aperiodic_rdegree),
                    std::max(periodic_ldegree, periodic_rdegree)>;
  using PeriodicPolynomial = typename Result::PeriodicPolynomial;

  auto aperiodic = left.aperiodic_ + right.aperiodic_;
  typename Result::PolynomialsByAngularFrequency periodic;
  periodic.reserve(left.periodic_.size() + right.periodic_.size());
  auto it_left = left.periodic_.cbegin();
  auto it_right = right.periodic_.cbegin();
  while (it_left != left.periodic_.cend() ||
         it_right != right.periodic_.cend()) {
    auto const ωl = it_left == left.periodic_.cend()
                        ? Infinity<AngularFrequency>
                        : it_left->first;
    auto const ωr = it_right == right.periodic_.cend()
                        ? Infinity<AngularFrequency>
                        : it_right->first;
    if (ωl < ωr) {
      auto const& polynomials_left = it_left->second;
      periodic.emplace_back(
          ωl,
          typename Result::Polynomials{
              .sin = PeriodicPolynomial(polynomials_left.sin),
              .cos = PeriodicPolynomial(polynomials_left.cos)});
      ++it_left;
    } else if (ωr < ωl) {
      auto const& polynomials_right = it_right->second;
      periodic.emplace_back(
          ωr,
          typename Result::Polynomials{
              .sin = PeriodicPolynomial(polynomials_right.sin),
              .cos = PeriodicPolynomial(polynomials_right.cos)});
      ++it_right;
    } else {
      DCHECK_EQ(ωl, ωr);
      auto const& polynomials_left = it_left->second;
      auto const& polynomials_right = it_right->second;
      periodic.emplace_back(
          ωl,
          typename Result::Polynomials{
              .sin = polynomials_left.sin + polynomials_right.sin,
              .cos = polynomials_left.cos + polynomials_right.cos});
      ++it_left;
      ++it_right;
    }
  }
  // Because we have done a merge on the periodic vectors, we can use the
  // trusted constructor.
  return {typename Result::TrustedPrivateConstructor{},
          std::move(aperiodic),
          std::move(periodic)};
}

template<typename Value,
         int aperiodic_ldegree, int periodic_ldegree,
         int aperiodic_rdegree, int periodic_rdegree>
PoissonSeries<Value,
              std::max(aperiodic_ldegree, aperiodic_rdegree),
              std::max(periodic_ldegree, periodic_rdegree)>
operator-(PoissonSeries<Value,
                        aperiodic_ldegree, periodic_ldegree> const& left,
          PoissonSeries<Value,
                        aperiodic_rdegree, periodic_rdegree> const& right) {
  using Result = PoissonSeries<Value,
                               std::max(aperiodic_ldegree, aperiodic_rdegree),
                               std::max(periodic_ldegree, periodic_rdegree)>;
  using PeriodicPolynomial = typename Result::PeriodicPolynomial;

  auto aperiodic = left.aperiodic_ - right.aperiodic_;
  typename Result::PolynomialsByAngularFrequency periodic;
  periodic.reserve(left.periodic_.size() + right.periodic_.size());
  auto it_left = left.periodic_.cbegin();
  auto it_right = right.periodic_.cbegin();
  while (it_left != left.periodic_.cend() ||
         it_right != right.periodic_.cend()) {
    auto const ωl = it_left == left.periodic_.cend()
                        ? Infinity<AngularFrequency>
                        : it_left->first;
    auto const ωr = it_right == right.periodic_.cend()
                        ? Infinity<AngularFrequency>
                        : it_right->first;
    if (ωl < ωr) {
      auto const& polynomials_left = it_left->second;
      periodic.emplace_back(
          ωl,
          typename Result::Polynomials{
              .sin = PeriodicPolynomial(polynomials_left.sin),
              .cos = PeriodicPolynomial(polynomials_left.cos)});
      ++it_left;
    } else if (ωr < ωl) {
      auto const& polynomials_right = it_right->second;
      periodic.emplace_back(
          ωr,
          typename Result::Polynomials{
              .sin = PeriodicPolynomial(-polynomials_right.sin),
              .cos = PeriodicPolynomial(-polynomials_right.cos)});
      ++it_right;
    } else {
      DCHECK_EQ(ωl, ωr);
      auto const& polynomials_left = it_left->second;
      auto const& polynomials_right = it_right->second;
      periodic.emplace_back(
          ωl,
          typename Result::Polynomials{
              .sin = polynomials_left.sin - polynomials_right.sin,
              .cos = polynomials_left.cos - polynomials_right.cos});
      ++it_left;
      ++it_right;
    }
  }
  // Because we have done a merge on the periodic vectors, we can use the
  // trusted constructor.
  return {typename Result::TrustedPrivateConstructor{},
          std::move(aperiodic),
          std::move(periodic)};
}

template<typename Scalar, typename Value,
         int aperiodic_rdegree, int periodic_rdegree>
PoissonSeries<Product<Scalar, Value>,
              aperiodic_rdegree, periodic_rdegree>
operator*(Scalar const& left,
          PoissonSeries<Value,
                        aperiodic_rdegree, periodic_rdegree> const& right) {
  using Result = PoissonSeries<Product<Scalar, Value>,
                               aperiodic_rdegree, periodic_rdegree>;
  auto aperiodic = left * right.aperiodic_;
  typename Result::PolynomialsByAngularFrequency periodic;
  periodic.reserve(right.periodic_.size());
  for (auto const& [ω, polynomials] : right.periodic_) {
    periodic.emplace_back(
        ω,
        typename Result::Polynomials{.sin = left * polynomials.sin,
                                     .cos = left * polynomials.cos});
  }
  return {typename Result::PrivateConstructor{},
          std::move(aperiodic),
          std::move(periodic)};
}

template<typename Scalar, typename Value,
         int aperiodic_ldegree, int periodic_ldegree>
PoissonSeries<Product<Value, Scalar>,
              aperiodic_ldegree, periodic_ldegree>
operator*(PoissonSeries<Value,
                        aperiodic_ldegree, periodic_ldegree> const& left,
          Scalar const& right) {
  using Result = PoissonSeries<Product<Value, Scalar>,
                               aperiodic_ldegree, periodic_ldegree>;
  auto aperiodic = left.aperiodic_ * right;
  typename Result::PolynomialsByAngularFrequency periodic;
  periodic.reserve(left.periodic_.size());
  for (auto const& [ω, polynomials] : left.periodic_) {
    periodic.emplace_back(
        ω,
        typename Result::Polynomials{.sin = polynomials.sin * right,
                                     .cos = polynomials.cos * right});
  }
  return {typename Result::TrustedPrivateConstructor{},
          std::move(aperiodic),
          std::move(periodic)};
}

template<typename Scalar, typename Value,
         int aperiodic_ldegree, int periodic_ldegree>
PoissonSeries<Quotient<Value, Scalar>,
              aperiodic_ldegree, periodic_ldegree>
operator/(PoissonSeries<Value,
                        aperiodic_ldegree, periodic_ldegree> const& left,
          Scalar const& right) {
  using Result = PoissonSeries<Quotient<Value, Scalar>,
                               aperiodic_ldegree, periodic_ldegree>;
  auto aperiodic = left.aperiodic_ / right;
  typename Result::PolynomialsByAngularFrequency periodic;
  periodic.reserve(left.periodic_.size());
  for (auto const& [ω, polynomials] : left.periodic_) {
    periodic.emplace_back(
        ω,
        typename Result::Polynomials{.sin = polynomials.sin / right,
                                     .cos = polynomials.cos / right});
  }
  return {typename Result::TrustedPrivateConstructor{},
          std::move(aperiodic),
          std::move(periodic)};
}

template<typename LValue, typename RValue,
         int aperiodic_ldegree, int periodic_ldegree,
         int aperiodic_rdegree, int periodic_rdegree>
PoissonSeries<Product<LValue, RValue>,
              std::max({aperiodic_ldegree + aperiodic_rdegree,
                        aperiodic_ldegree + periodic_rdegree,
                        periodic_ldegree + aperiodic_rdegree,
                        periodic_ldegree + periodic_rdegree}),
              std::max({aperiodic_ldegree + periodic_rdegree,
                        periodic_ldegree + aperiodic_rdegree,
                        periodic_ldegree + periodic_rdegree})>
operator*(PoissonSeries<LValue,
                        aperiodic_ldegree, periodic_ldegree> const& left,
          PoissonSeries<RValue,
                        aperiodic_rdegree, periodic_rdegree> const& right) {
  auto product = [](auto const& left, auto const& right) {
    return left * right;
  };

  return Multiply<LValue, RValue,
                  aperiodic_ldegree, periodic_ldegree,
                  aperiodic_rdegree, periodic_rdegree>(left, right, product);
}

template<typename LValue, typename RValue,
         int aperiodic_ldegree, int periodic_ldegree,
         int aperiodic_rdegree, int periodic_rdegree>
PoissonSeries<typename Hilbert<LValue, RValue>::InnerProductType,
              std::max({aperiodic_ldegree + aperiodic_rdegree,
                        aperiodic_ldegree + periodic_rdegree,
                        periodic_ldegree + aperiodic_rdegree,
                        periodic_ldegree + periodic_rdegree}),
              std::max({aperiodic_ldegree + periodic_rdegree,
                        periodic_ldegree + aperiodic_rdegree,
                        periodic_ldegree + periodic_rdegree})>
PointwiseInnerProduct(
    PoissonSeries<LValue,
                  aperiodic_ldegree, periodic_ldegree> const& left,
    PoissonSeries<RValue,
                  aperiodic_rdegree, periodic_rdegree> const& right) {
  auto product = [](auto const& left, auto const& right) {
    return PointwiseInnerProduct(left, right);
  };

  return Multiply<LValue, RValue,
                  aperiodic_ldegree,  periodic_ldegree,
                  aperiodic_rdegree, periodic_rdegree>(left, right, product);
}

template<typename Value,
         int aperiodic_degree, int periodic_degree>
std::ostream& operator<<(
    std::ostream& out,
    PoissonSeries<Value,
                  aperiodic_degree, periodic_degree> const& series) {
  bool is_start_of_output = true;
  if (!series.aperiodic_.is_zero()) {
    out << series.aperiodic_;
    is_start_of_output = false;
  }
  for (auto const& [ω, polynomials] : series.periodic_) {
    if (!polynomials.sin.is_zero()) {
      if (!is_start_of_output) {
        out << " + ";
      }
      out <<"(" << polynomials.sin << ") * Sin(" << DebugString(ω)
          << " * (T - " << series.origin_ << "))";
      is_start_of_output = false;
    }
    if (!polynomials.cos.is_zero()) {
      if (!is_start_of_output) {
        out << " + ";
      }
      out << "(" << polynomials.cos << ") * Cos(" << DebugString(ω)
          << " * (T - " << series.origin_ << "))";
      is_start_of_output = false;
    }
  }
  return out;
}

template<typename LValue, typename RValue,
         int aperiodic_ldegree, int periodic_ldegree,
         int aperiodic_rdegree, int periodic_rdegree,
         int aperiodic_wdegree, int periodic_wdegree>
typename Hilbert<LValue, RValue>::InnerProductType InnerProduct(
    PoissonSeries<LValue,
                  aperiodic_ldegree, periodic_ldegree> const& left,
    PoissonSeries<RValue,
                  aperiodic_rdegree, periodic_rdegree> const& right,
    PoissonSeries<double,
                  aperiodic_wdegree, periodic_wdegree> const& weight,
    Instant const& t_min,
    Instant const& t_max) {
  AngularFrequency const ω_cutoff =
      2 * π * Radian * clenshaw_curtis_max_periods_overall / (t_max - t_min);
  auto const left_split = left.Split(ω_cutoff);
  auto const right_split = right.Split(ω_cutoff);

  AngularFrequency const max_ω =
      left_split.slow.max_ω() + right_split.slow.max_ω() + weight.max_ω();
  std::optional<int> const max_points =
      _quadrature::MaxPointsHeuristicsForAutomaticClenshawCurtis(
          max_ω,
          t_max - t_min,
          clenshaw_curtis_min_points_overall,
          clenshaw_curtis_points_per_period);

  auto slow_integrand = [&left_split, &right_split, &weight](Instant const& t) {
    return Hilbert<LValue, RValue>::InnerProduct(left_split.slow(t),
                                                 right_split.slow(t)) *
           weight(t);
  };
  auto const slow_quadrature = _quadrature::AutomaticClenshawCurtis(
      slow_integrand,
      t_min, t_max,
      /*max_relative_error=*/clenshaw_curtis_relative_error,
      /*max_points=*/max_points);

  auto const fast_integrand =
      (PointwiseInnerProduct(left_split.fast, right_split.slow) +
       PointwiseInnerProduct(left_split.slow, right_split.fast) +
       PointwiseInnerProduct(left_split.fast, right_split.fast)) *
      weight;
  auto const fast_quadrature = fast_integrand.Integrate(t_min, t_max);

  return (slow_quadrature + fast_quadrature) / (t_max - t_min);
}

}  // namespace internal
}  // namespace _poisson_series
}  // namespace numerics
}  // namespace principia
