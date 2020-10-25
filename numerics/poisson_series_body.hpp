
#pragma once

#include "numerics/poisson_series.hpp"

#include <algorithm>
#include <functional>
#include <map>
#include <optional>
#include <utility>
#include <vector>

#include "numerics/double_precision.hpp"
#include "numerics/quadrature.hpp"
#include "numerics/ulp_distance.hpp"
#include "quantities/elementary_functions.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace numerics {
namespace internal_poisson_series {

using quantities::Abs;
using quantities::Cos;
using quantities::Infinity;
using quantities::Sin;
using quantities::Sqrt;
using quantities::Variation;
using quantities::si::Radian;
namespace si = quantities::si;

// These parameters have been tuned for approximation of the Moon over 3 months
// with 10 periods.

// In the fast/slow algorithms, specifies the maximum number of periods over the
// time interval below which we use Clenshaw-Curtis integration.
constexpr int clenshaw_curtis_max_periods_overall = 40;

// The minimum value of the max_point parameter passed to Clenshaw-Curtis
// integration, irrespective of the frequencies of the argument function.
constexpr int clenshaw_curtis_min_points_overall = 33;

// The maximum number of points use in Clenshaw-Curtis integration for each
// period of the highest frequency of the argument function.
constexpr int clenshaw_curtis_point_per_period = 4;

// The desired relative error on Clenshaw-Curtis integration, as determined by
// two successive computations with increasing number of points.
constexpr double clenshaw_curtis_relative_error = 0x1p-32;

template<typename Value, int degree_,
         template<typename, typename, int> class Evaluator>
typename PoissonSeries<Primitive<Value, Time>, degree_ + 1, Evaluator>::
    Polynomials
AngularFrequencyPrimitive(
    AngularFrequency const& ω,
    typename PoissonSeries<Value, degree_, Evaluator>::Polynomials const&
        polynomials) {
  using Result = PoissonSeries<Primitive<Value, Time>, degree_ + 1, Evaluator>;

  // Integration by parts.
  typename Result::Polynomials const first_part{
      /*sin=*/typename Result::Polynomial(polynomials.cos / ω * Radian),
      /*cos=*/typename Result::Polynomial(-polynomials.sin / ω * Radian)};
  if constexpr (degree_ == 0) {
    return first_part;
  } else {
    auto const sin_polynomial =
        -polynomials.cos.template Derivative<1>() / ω * Radian;
    auto const cos_polynomial =
        polynomials.sin.template Derivative<1>() / ω * Radian;
    auto const second_part =
        AngularFrequencyPrimitive<Value, degree_ - 1, Evaluator>(
            ω,
            {/*sin=*/sin_polynomial,
             /*cos=*/cos_polynomial});
    return {/*sin=*/first_part.sin + second_part.sin,
            /*cos=*/first_part.cos + second_part.cos};
  }
}

// A helper for multiplication of Poisson series and pointwise inner product.
// The functor Product must take a pair of Poisson series with the types of left
// and right and return a suitable Poisson series.
template<typename LValue, typename RValue,
         int ldegree_, int rdegree_,
         template<typename, typename, int> class Evaluator,
         typename Product>
auto Multiply(PoissonSeries<LValue, ldegree_, Evaluator> const& left,
              PoissonSeries<RValue, rdegree_, Evaluator> const& right,
              Product const& product) {
  using Result = PoissonSeries<
      typename std::invoke_result_t<
          Product,
          typename PoissonSeries<LValue, ldegree_, Evaluator>::Polynomial,
          typename PoissonSeries<RValue, rdegree_, Evaluator>::Polynomial>::
          Value,
      ldegree_ + rdegree_,
      Evaluator>;

  auto aperiodic = product(left.aperiodic_, right.aperiodic_);

  // Compute all the individual terms using elementary trigonometric identities
  // and put them in a vector, because the same frequency may appear multiple
  // times.
  typename Result::PolynomialsByAngularFrequency periodic;
  periodic.reserve(left.periodic_.size() + right.periodic_.size() +
                   2 * left.periodic_.size() * right.periodic_.size());
  for (auto const& [ω, polynomials] : left.periodic_) {
    periodic.emplace_back(
        ω,
        typename Result::Polynomials{
            /*sin=*/product(polynomials.sin, right.aperiodic_),
            /*cos=*/product(polynomials.cos, right.aperiodic_)});
  }
  for (auto const& [ω, polynomials] : right.periodic_) {
    periodic.emplace_back(
        ω,
        typename Result::Polynomials{
            /*sin=*/product(left.aperiodic_, polynomials.sin),
            /*cos=*/product(left.aperiodic_, polynomials.cos)});
  }
  for (auto const& [ωl, polynomials_left] : left.periodic_) {
    for (auto const& [ωr, polynomials_right] : right.periodic_) {
      auto const cos_cos = product(polynomials_left.cos, polynomials_right.cos);
      auto const cos_sin = product(polynomials_left.cos, polynomials_right.sin);
      auto const sin_cos = product(polynomials_left.sin, polynomials_right.cos);
      auto const sin_sin = product(polynomials_left.sin, polynomials_right.sin);
      periodic.emplace_back(
          ωl - ωr,
          typename Result::Polynomials{/*sin=*/(-cos_sin + sin_cos) / 2,
                                       /*cos=*/(sin_sin + cos_cos) / 2});
      periodic.emplace_back(
          ωl + ωr,
          typename Result::Polynomials{/*sin=*/(cos_sin + sin_cos) / 2,
                                       /*cos=*/(-sin_sin + cos_cos) / 2});
    }
  }

  return Result(typename Result::PrivateConstructor{},
                std::move(aperiodic),
                std::move(periodic));
}

template<typename Value, int degree_,
         template<typename, typename, int> class Evaluator>
PoissonSeries<Value, degree_, Evaluator>::PoissonSeries(
    Polynomial const& aperiodic,
    PolynomialsByAngularFrequency const& periodic)
    : PoissonSeries(PrivateConstructor{},
                    Polynomial(aperiodic),
                    PolynomialsByAngularFrequency(periodic)) {}

template<typename Value, int degree_,
         template<typename, typename, int> class Evaluator>
template<int higher_degree_,
         template<typename, typename, int> class HigherEvaluator>
PoissonSeries<Value, degree_, Evaluator>::
operator PoissonSeries<Value, higher_degree_, HigherEvaluator>() const {
  static_assert(degree_ <= higher_degree_);
  using Result = PoissonSeries<Value, higher_degree_, HigherEvaluator>;
  auto aperiodic = typename Result::Polynomial(aperiodic_);
  typename Result::PolynomialsByAngularFrequency periodic;
  periodic.reserve(periodic_.size());
  for (auto const& [ω, polynomials] : periodic_) {
    periodic.emplace_back(
        ω,
        typename Result::Polynomials{
            /*sin=*/typename Result::Polynomial(polynomials.sin),
            /*cos=*/typename Result::Polynomial(polynomials.cos)});
  }
  return Result(typename Result::TrustedPrivateConstructor{},
                std::move(aperiodic),
                std::move(periodic));
}

template<typename Value, int degree_,
         template<typename, typename, int> class Evaluator>
Instant const& PoissonSeries<Value, degree_, Evaluator>::origin() const {
  return origin_;
}

template<typename Value, int degree_,
         template<typename, typename, int> class Evaluator>
Value PoissonSeries<Value, degree_, Evaluator>::operator()(
    Instant const& t) const {
  Value result = aperiodic_(t);
  for (auto const& [ω, polynomials] : periodic_) {
    result += polynomials.sin(t) * Sin(ω * (t - origin_)) +
              polynomials.cos(t) * Cos(ω * (t - origin_));
  }
  return result;
}

template<typename Value, int degree_,
         template<typename, typename, int> class Evaluator>
PoissonSeries<Value, degree_, Evaluator>
PoissonSeries<Value, degree_, Evaluator>::AtOrigin(
    Instant const& origin) const {
  Time const shift = origin - origin_;
  auto aperiodic = aperiodic_.AtOrigin(origin);

  PolynomialsByAngularFrequency periodic;
  periodic.reserve(periodic_.size());
  for (auto const& [ω, polynomials] : periodic_) {
    double const cos_ω_shift = Cos(ω * shift);
    double const sin_ω_shift = Sin(ω * shift);
    Polynomial const sin_at_origin = polynomials.sin.AtOrigin(origin);
    Polynomial const cos_at_origin = polynomials.cos.AtOrigin(origin);
    periodic.emplace_back(ω,
                          Polynomials{/*sin=*/sin_at_origin * cos_ω_shift -
                                              cos_at_origin * sin_ω_shift,
                                      /*cos=*/sin_at_origin * sin_ω_shift +
                                              cos_at_origin * cos_ω_shift});
  }
  return {TrustedPrivateConstructor{},
          std::move(aperiodic),
          std::move(periodic)};
}

template<typename Value, int degree_,
         template<typename, typename, int> class Evaluator>
PoissonSeries<quantities::Primitive<Value, Time>, degree_ + 1, Evaluator>
PoissonSeries<Value, degree_, Evaluator>::Primitive() const {
  using Result =
      PoissonSeries<quantities::Primitive<Value, Time>, degree_ + 1, Evaluator>;
  typename Result::Polynomial aperiodic = aperiodic_.Primitive();
  typename Result::PolynomialsByAngularFrequency periodic;
  periodic.reserve(periodic_.size());
  for (auto const& [ω, polynomials] : periodic_) {
    periodic.emplace_back(
        ω,
        AngularFrequencyPrimitive<Value, degree_, Evaluator>(ω, polynomials));
  }
  return Result{aperiodic, periodic};
}

template<typename Value, int degree_,
         template<typename, typename, int> class Evaluator>
quantities::Primitive<Value, Time>
PoissonSeries<Value, degree_, Evaluator>::Integrate(Instant const& t1,
                                                    Instant const& t2) const {
  using FirstPart = PoissonSeries<Value, degree_, Evaluator>;
  using SecondPart = PoissonSeries<Variation<Value>, degree_ - 1, Evaluator>;
  auto const aperiodic_primitive = aperiodic_.Primitive();
  quantities::Primitive<Value, Time> result =
      aperiodic_primitive(t2) - aperiodic_primitive(t1);
  for (auto const& [ω, polynomials] : periodic_) {
    // Integration by parts.
    FirstPart const first_part(
        typename FirstPart::Polynomial({}, origin_),
        {{ω,
          {/*sin=*/typename FirstPart::Polynomial(polynomials.cos),
           /*cos=*/typename FirstPart::Polynomial(-polynomials.sin)}}});
    DoublePrecision<Value> sum;
    sum += first_part(t2);
    sum -= first_part(t1);

    if constexpr (degree_ != 0) {
      auto const sin_polynomial =
          -polynomials.cos.template Derivative<1>();
      auto const cos_polynomial =
          polynomials.sin.template Derivative<1>();
      SecondPart const second_part(typename SecondPart::Polynomial({}, origin_),
                                   {{ω,
                                     {/*sin=*/sin_polynomial,
                                      /*cos=*/cos_polynomial}}});
      sum += second_part.Integrate(t1, t2);
    }
    result += (sum.value + sum.error) / ω * Radian;
  }
  return result;
}

template<typename Value, int degree_,
         template<typename, typename, int> class Evaluator>
template<int wdegree_>
typename Hilbert<Value>::NormType
PoissonSeries<Value, degree_, Evaluator>::Norm(
    PoissonSeries<double, wdegree_, Evaluator> const& weight,
    Instant const& t_min,
    Instant const& t_max) const {
  AngularFrequency const ω_cutoff =
      2 * π * Radian * clenshaw_curtis_max_periods_overall / (t_max - t_min);
  auto const split = Split(ω_cutoff);

  AngularFrequency const max_ω =
      (split.slow.periodic_.empty() ? AngularFrequency{}
                                     : 2 * split.slow.periodic_.back().first) +
      (weight.periodic_.empty() ? AngularFrequency{}
                                : weight.periodic_.back().first);
  std::optional<int> max_points =
      max_ω == AngularFrequency()
          ? std::optional<int>{}
          : std::max(
                clenshaw_curtis_min_points_overall,
                static_cast<int>(clenshaw_curtis_point_per_period *
                                 (t_max - t_min) * max_ω / (2 * π * Radian)));

  auto slow_integrand = [&split, &weight](Instant const& t) {
    return Hilbert<Value>::Norm²(split.slow(t)) * weight(t);
  };
  auto const slow_quadrature = quadrature::AutomaticClenshawCurtis(
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

template<typename Value, int degree_,
         template<typename, typename, int> class Evaluator>
template<int d>
PoissonSeries<Value, degree_, Evaluator>&
PoissonSeries<Value, degree_, Evaluator>::operator+=(
    PoissonSeries<Value, d, Evaluator> const& right) {
  static_assert(d <= degree);
  *this = *this + right;
  return *this;
}

template<typename Value, int degree_,
         template<typename, typename, int> class Evaluator>
template<int d>
PoissonSeries<Value, degree_, Evaluator>&
PoissonSeries<Value, degree_, Evaluator>::operator-=(
    PoissonSeries<Value, d, Evaluator> const& right) {
  static_assert(d <= degree);
  *this = *this - right;
  return *this;
}

template<typename Value, int degree_,
         template<typename, typename, int> class Evaluator>
void PoissonSeries<Value, degree_, Evaluator>::WriteToMessage(
    not_null<serialization::PoissonSeries*> const message) const {
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

template<typename Value, int degree_,
         template<typename, typename, int> class Evaluator>
PoissonSeries<Value, degree_, Evaluator>
PoissonSeries<Value, degree_, Evaluator>::ReadFromMessage(
    serialization::PoissonSeries const& message) {
  auto const aperiodic = Polynomial::ReadFromMessage(message.aperiodic());
  PolynomialsByAngularFrequency periodic;
  for (auto const& polynomial_and_angular_frequency : message.periodic()) {
    auto const ω = AngularFrequency::ReadFromMessage(
        polynomial_and_angular_frequency.angular_frequency());
    auto const sin =
        Polynomial::ReadFromMessage(polynomial_and_angular_frequency.sin());
    auto const cos =
        Polynomial::ReadFromMessage(polynomial_and_angular_frequency.cos());
    periodic.emplace_back(ω, Polynomials{/*sin=*/sin, /*cos=*/cos});
  }
  return PoissonSeries(aperiodic, periodic);
}

template<typename Value, int degree_,
         template<typename, typename, int> class Evaluator>
PoissonSeries<Value, degree_, Evaluator>::PoissonSeries(
    PrivateConstructor,
    Polynomial aperiodic,
    PolynomialsByAngularFrequency periodic)
    : origin_(aperiodic.origin()),
      aperiodic_(std::move(aperiodic)),
      periodic_(std::move(periodic)) {
  // The |periodic| vector may have elements with positive or negative angular
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
    } else {
      aperiodic_ += polynomials.cos;
      it = periodic_.erase(it);
    }
    previous_abs_ω = abs_ω;
  }
}

template<typename Value, int degree_,
         template<typename, typename, int> class Evaluator>
PoissonSeries<Value, degree_, Evaluator>::PoissonSeries(
    TrustedPrivateConstructor,
    Polynomial aperiodic,
    PolynomialsByAngularFrequency periodic)
    : origin_(aperiodic.origin()),
      aperiodic_(std::move(aperiodic)),
      periodic_(std::move(periodic)) {}

template<typename Value, int degree_,
         template<typename, typename, int> class Evaluator>
typename PoissonSeries<Value, degree_, Evaluator>::SplitPoissonSeries
PoissonSeries<Value, degree_, Evaluator>::Split(
    AngularFrequency const& ω_cutoff) const {
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
                       typename PoissonSeries::Polynomial(
                           typename PoissonSeries::Polynomial::Coefficients{},
                           aperiodic_.origin()),
                       std::move(slow_periodic));
    PoissonSeries fast(TrustedPrivateConstructor{},
                       aperiodic_,
                       std::move(fast_periodic));
    return {/*slow=*/std::move(slow), /*fast=*/std::move(fast)};
  } else {
    PoissonSeries slow(TrustedPrivateConstructor{},
                       aperiodic_,
                       std::move(slow_periodic));
    PoissonSeries fast(TrustedPrivateConstructor{},
                       typename PoissonSeries::Polynomial(
                           typename PoissonSeries::Polynomial::Coefficients{},
                           aperiodic_.origin()),
                       std::move(fast_periodic));
    return {/*slow=*/std::move(slow), /*fast=*/std::move(fast)};
  }
}

template<typename Value, int rdegree_,
         template<typename, typename, int> class Evaluator>
PoissonSeries<Value, rdegree_, Evaluator> operator+(
    PoissonSeries<Value, rdegree_, Evaluator> const& right) {
  return right;
}

template<typename Value, int rdegree_,
         template<typename, typename, int> class Evaluator>
PoissonSeries<Value, rdegree_, Evaluator>
operator-(PoissonSeries<Value, rdegree_, Evaluator> const& right) {
  using Result = PoissonSeries<Value, rdegree_, Evaluator>;
  auto aperiodic = -right.aperiodic_;
  typename Result::PolynomialsByAngularFrequency periodic;
  periodic.reserve(right.periodic_.size());
  for (auto const& [ω, polynomials] : right.periodic_) {
    periodic.emplace_back(
        ω,
        typename Result::Polynomials{/*sin=*/-polynomials.sin,
                                     /*cos=*/-polynomials.cos});
  }
  return {typename Result::TrustedPrivateConstructor{},
          std::move(aperiodic),
          std::move(periodic)};
}

template<typename Value, int ldegree_, int rdegree_,
         template<typename, typename, int> class Evaluator>
PoissonSeries<Value, std::max(ldegree_, rdegree_), Evaluator>
operator+(PoissonSeries<Value, ldegree_, Evaluator> const& left,
          PoissonSeries<Value, rdegree_, Evaluator> const& right) {
  using Result = PoissonSeries<Value, std::max(ldegree_, rdegree_), Evaluator>;
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
              /*sin=*/typename Result::Polynomial(polynomials_left.sin),
              /*cos=*/typename Result::Polynomial(polynomials_left.cos)});
      ++it_left;
    } else if (ωr < ωl) {
      auto const& polynomials_right = it_right->second;
      periodic.emplace_back(
          ωr,
          typename Result::Polynomials{
              /*sin=*/typename Result::Polynomial(polynomials_right.sin),
              /*cos=*/typename Result::Polynomial(polynomials_right.cos)});
      ++it_right;
    } else {
      DCHECK_EQ(ωl, ωr);
      auto const& polynomials_left = it_left->second;
      auto const& polynomials_right = it_right->second;
      periodic.emplace_back(
          ωl,
          typename Result::Polynomials{
              /*sin=*/polynomials_left.sin + polynomials_right.sin,
              /*cos=*/polynomials_left.cos + polynomials_right.cos});
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

template<typename Value, int ldegree_, int rdegree_,
         template<typename, typename, int> class Evaluator>
PoissonSeries<Value, std::max(ldegree_, rdegree_), Evaluator>
operator-(PoissonSeries<Value, ldegree_, Evaluator> const& left,
          PoissonSeries<Value, rdegree_, Evaluator> const& right) {
  using Result = PoissonSeries<Value, std::max(ldegree_, rdegree_), Evaluator>;
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
              /*sin=*/typename Result::Polynomial(polynomials_left.sin),
              /*cos=*/typename Result::Polynomial(polynomials_left.cos)});
      ++it_left;
    } else if (ωr < ωl) {
      auto const& polynomials_right = it_right->second;
      periodic.emplace_back(
          ωr,
          typename Result::Polynomials{
              /*sin=*/typename Result::Polynomial(-polynomials_right.sin),
              /*cos=*/typename Result::Polynomial(-polynomials_right.cos)});
      ++it_right;
    } else {
      DCHECK_EQ(ωl, ωr);
      auto const& polynomials_left = it_left->second;
      auto const& polynomials_right = it_right->second;
      periodic.emplace_back(
          ωl,
          typename Result::Polynomials{
              /*sin=*/polynomials_left.sin - polynomials_right.sin,
              /*cos=*/polynomials_left.cos - polynomials_right.cos});
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

template<typename Scalar, typename Value, int degree_,
         template<typename, typename, int> class Evaluator>
PoissonSeries<Product<Scalar, Value>, degree_, Evaluator>
operator*(Scalar const& left,
          PoissonSeries<Value, degree_, Evaluator> const& right) {
  using Result = PoissonSeries<Product<Scalar, Value>, degree_, Evaluator>;
  auto aperiodic = left * right.aperiodic_;
  typename Result::PolynomialsByAngularFrequency periodic;
  periodic.reserve(right.periodic_.size());
  for (auto const& [ω, polynomials] : right.periodic_) {
    periodic.emplace_back(
        ω,
        typename Result::Polynomials{/*sin=*/left * polynomials.sin,
                                     /*cos=*/left * polynomials.cos});
  }
  return {typename Result::PrivateConstructor{},
          std::move(aperiodic),
          std::move(periodic)};
}

template<typename Scalar, typename Value, int degree_,
         template<typename, typename, int> class Evaluator>
PoissonSeries<Product<Value, Scalar>, degree_, Evaluator>
operator*(PoissonSeries<Value, degree_, Evaluator> const& left,
          Scalar const& right) {
  using Result = PoissonSeries<Product<Value, Scalar>, degree_, Evaluator>;
  auto aperiodic = left.aperiodic_ * right;
  typename Result::PolynomialsByAngularFrequency periodic;
  periodic.reserve(left.periodic_.size());
  for (auto const& [ω, polynomials] : left.periodic_) {
    periodic.emplace_back(
        ω,
        typename Result::Polynomials{/*sin=*/polynomials.sin * right,
                                     /*cos=*/polynomials.cos * right});
  }
  return {typename Result::TrustedPrivateConstructor{},
          std::move(aperiodic),
          std::move(periodic)};
}

template<typename Scalar, typename Value, int degree_,
         template<typename, typename, int> class Evaluator>
PoissonSeries<Quotient<Value, Scalar>, degree_, Evaluator>
operator/(PoissonSeries<Value, degree_, Evaluator> const& left,
          Scalar const& right) {
  using Result = PoissonSeries<Quotient<Value, Scalar>, degree_, Evaluator>;
  auto aperiodic = left.aperiodic_ / right;
  typename Result::PolynomialsByAngularFrequency periodic;
  periodic.reserve(left.periodic_.size());
  for (auto const& [ω, polynomials] : left.periodic_) {
    periodic.emplace_back(
        ω,
        typename Result::Polynomials{/*sin=*/polynomials.sin / right,
                                     /*cos=*/polynomials.cos / right});
  }
  return {typename Result::TrustedPrivateConstructor{},
          std::move(aperiodic),
          std::move(periodic)};
}

template<typename LValue, typename RValue,
         int ldegree_, int rdegree_,
         template<typename, typename, int> class Evaluator>
PoissonSeries<Product<LValue, RValue>, ldegree_ + rdegree_, Evaluator>
operator*(PoissonSeries<LValue, ldegree_, Evaluator> const& left,
          PoissonSeries<RValue, rdegree_, Evaluator> const& right) {
  auto product =
      [](typename PoissonSeries<LValue, ldegree_, Evaluator>::Polynomial const&
             left,
         typename PoissonSeries<RValue, rdegree_, Evaluator>::Polynomial const&
             right) {
    return left * right;
  };

  return Multiply<LValue, RValue, ldegree_, rdegree_, Evaluator>(
      left, right, product);
}

template<typename LValue, typename RValue,
         int ldegree_, int rdegree_,
         template<typename, typename, int> class Evaluator>
PoissonSeries<typename Hilbert<LValue, RValue>::InnerProductType,
              ldegree_ + rdegree_,
              Evaluator>
PointwiseInnerProduct(PoissonSeries<LValue, ldegree_, Evaluator> const& left,
                      PoissonSeries<RValue, rdegree_, Evaluator> const& right) {
  auto product =
      [](typename PoissonSeries<LValue, ldegree_, Evaluator>::Polynomial const&
             left,
         typename PoissonSeries<RValue, rdegree_, Evaluator>::Polynomial const&
             right) {
    return PointwiseInnerProduct(left, right);
  };

  return Multiply<LValue, RValue, ldegree_, rdegree_, Evaluator>(
      left, right, product);
}

template<typename Value, int degree_,
         template<typename, typename, int> class Evaluator>
std::ostream& operator<<(
    std::ostream& out,
    PoissonSeries<Value, degree_, Evaluator> const& series) {
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
      out <<"(" << polynomials.sin << ") * Sin(" << quantities::DebugString(ω)
          << " * (T - " << series.origin_ << "))";
      is_start_of_output = false;
    }
    if (!polynomials.cos.is_zero()) {
      if (!is_start_of_output) {
        out << " + ";
      }
      out << "(" << polynomials.cos << ") * Cos(" << quantities::DebugString(ω)
          << " * (T - " << series.origin_ << "))";
      is_start_of_output = false;
    }
  }
  return out;
}

template<typename LValue, typename RValue,
         int ldegree_, int rdegree_, int wdegree_,
         template<typename, typename, int> class Evaluator>
typename Hilbert<LValue, RValue>::InnerProductType InnerProduct(
    PoissonSeries<LValue, ldegree_, Evaluator> const& left,
    PoissonSeries<RValue, rdegree_, Evaluator> const& right,
    PoissonSeries<double, wdegree_, Evaluator> const& weight,
    Instant const& t_min,
    Instant const& t_max) {
  AngularFrequency const ω_cutoff =
      2 * π * Radian * clenshaw_curtis_max_periods_overall / (t_max - t_min);
  auto const left_split = left.Split(ω_cutoff);
  auto const right_split = right.Split(ω_cutoff);

  AngularFrequency const max_ω =
      (left_split.slow.periodic_.empty()
           ? AngularFrequency{}
           : left_split.slow.periodic_.back().first) +
      (right_split.slow.periodic_.empty()
           ? AngularFrequency{}
           : right_split.slow.periodic_.back().first) +
      (weight.periodic_.empty() ? AngularFrequency{}
                                : weight.periodic_.back().first);
  std::optional<int> max_points =
      max_ω == AngularFrequency()
          ? std::optional<int>{}
          : std::max(
                clenshaw_curtis_min_points_overall,
                static_cast<int>(clenshaw_curtis_point_per_period *
                                 (t_max - t_min) * max_ω / (2 * π * Radian)));

  auto slow_integrand = [&left_split, &right_split, &weight](Instant const& t) {
    return Hilbert<LValue, RValue>::InnerProduct(left_split.slow(t),
                                                 right_split.slow(t)) *
           weight(t);
  };
  auto const slow_quadrature = quadrature::AutomaticClenshawCurtis(
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

}  // namespace internal_poisson_series
}  // namespace numerics
}  // namespace principia
