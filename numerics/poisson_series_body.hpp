
#pragma once

#include "numerics/poisson_series.hpp"

#include <algorithm>
#include <functional>
#include <map>
#include <optional>

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
using quantities::Primitive;
using quantities::Sin;
using quantities::Variation;
using quantities::si::Radian;
using quantities::si::Second;
namespace si = quantities::si;

template<typename Value, int degree_,
         template<typename, typename, int> class Evaluator>
typename PoissonSeries<Primitive<Value, Time>, degree_ + 1, Evaluator>::
    Polynomials
AngularFrequencyPrimitive(
    AngularFrequency const& ω,
    typename PoissonSeries<Value, degree_, Evaluator>::Polynomials const&
        polynomials) {
  using Argument = PoissonSeries<Value, degree_, Evaluator>;
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

template<typename Value, int degree_,
         template<typename, typename, int> class Evaluator>
PoissonSeries<Value, degree_, Evaluator>::PoissonSeries(
    Polynomial const& aperiodic,
    PolynomialsByAngularFrequency const& periodic)
    : origin_(aperiodic.origin()),
      aperiodic_(aperiodic) {
  // The |periodic| map may have elements with positive or negative angular
  // frequencies.  Normalize our member variable to only have positive angular
  // frequencies.
  for (auto it = periodic.crbegin(); it != periodic.crend(); ++it) {
    auto const ω = it->first;
    auto const polynomials = it->second;

    // All polynomials must have the same origin.
    CHECK_EQ(origin_, polynomials.sin.origin());
    CHECK_EQ(origin_, polynomials.cos.origin());

    if (ω < AngularFrequency{}) {
      auto const positive_it = periodic_.find(-ω);
      if (positive_it == periodic_.cend()) {
        periodic_.emplace(-ω,
                          Polynomials{/*sin=*/-polynomials.sin,
                                      /*cos=*/polynomials.cos});
      } else {
        positive_it->second.sin -= polynomials.sin;
        positive_it->second.cos += polynomials.cos;
      }
    } else if (ω > AngularFrequency{}) {
      periodic_.insert(periodic_.cbegin(), *it);
    } else {
      aperiodic_ += polynomials.cos;
    }
  }
}

template<typename Value, int degree_,
         template<typename, typename, int> class Evaluator>
Value PoissonSeries<Value, degree_, Evaluator>::Evaluate(
    Instant const& t) const {
  Value result = aperiodic_.Evaluate(t);
  for (auto const& [ω, polynomials] : periodic_) {
    result += polynomials.sin.Evaluate(t) * Sin(ω * (t - origin_)) +
              polynomials.cos.Evaluate(t) * Cos(ω * (t - origin_));
  }
  return result;
}

template<typename Value, int degree_,
         template<typename, typename, int> class Evaluator>
PoissonSeries<quantities::Primitive<Value, Time>, degree_ + 1, Evaluator>
PoissonSeries<Value, degree_, Evaluator>::Primitive() const {
  using Result =
      PoissonSeries<quantities::Primitive<Value, Time>, degree_ + 1, Evaluator>;
  typename Result::PolynomialsByAngularFrequency periodic;
  typename Result::Polynomial const aperiodic = aperiodic_.Primitive();
  for (auto const& [ω, polynomials] : periodic_) {
    periodic.emplace_hint(
        periodic.cend(),
        ω,
        AngularFrequencyPrimitive<Value, degree_, Evaluator>(ω, polynomials));
  }
  return Result{aperiodic, std::move(periodic)};
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
  typename Result::PolynomialsByAngularFrequency periodic;
  auto const aperiodic = -right.aperiodic_;
  for (auto const& [ω, polynomials] : right.periodic_) {
    periodic.emplace_hint(
        periodic.cend(),
        ω,
        typename Result::Polynomials{/*sin=*/-polynomials.sin,
                                     /*cos=*/-polynomials.cos});
  }
  return {aperiodic, std::move(periodic)};
}

template<typename Value, int ldegree_, int rdegree_,
         template<typename, typename, int> class Evaluator>
PoissonSeries<Value, std::max(ldegree_, rdegree_), Evaluator>
operator+(PoissonSeries<Value, ldegree_, Evaluator> const& left,
          PoissonSeries<Value, rdegree_, Evaluator> const& right) {
  using Result = PoissonSeries<Value, std::max(ldegree_, rdegree_), Evaluator>;
  typename Result::PolynomialsByAngularFrequency periodic;
  auto const aperiodic = left.aperiodic_ + right.aperiodic_;
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
      periodic.insert(periodic.cend(), *it_left);
      ++it_left;
    } else if (ωr < ωl) {
      periodic.insert(periodic.cend(), *it_right);
      ++it_right;
    } else {
      DCHECK_EQ(ωl, ωr);
      auto const& polynomials_left = it_left->second;
      auto const& polynomials_right = it_right->second;
      periodic.emplace_hint(
          periodic.cend(),
          ωl,
          typename Result::Polynomials{
              /*sin=*/polynomials_left.sin + polynomials_right.sin,
              /*cos=*/polynomials_left.cos + polynomials_right.cos});
      ++it_left;
      ++it_right;
    }
  }
  return {aperiodic, std::move(periodic)};
}

template<typename Value, int ldegree_, int rdegree_,
         template<typename, typename, int> class Evaluator>
PoissonSeries<Value, std::max(ldegree_, rdegree_), Evaluator>
operator-(PoissonSeries<Value, ldegree_, Evaluator> const& left,
          PoissonSeries<Value, rdegree_, Evaluator> const& right) {
  using Result = PoissonSeries<Value, std::max(ldegree_, rdegree_), Evaluator>;
  typename Result::PolynomialsByAngularFrequency periodic;
  auto const aperiodic = left.aperiodic_ - right.aperiodic_;
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
      periodic.insert(periodic.cend(), *it_left);
      ++it_left;
    } else if (ωr < ωl) {
      auto const& polynomials_right = it_right->second;
      periodic.emplace_hint(
          periodic.cend(),
          ωr,
          typename Result::Polynomials{/*sin=*/-polynomials_right.sin,
                                       /*cos=*/-polynomials_right.cos});
      ++it_right;
    } else {
      DCHECK_EQ(ωl, ωr);
      auto const& polynomials_left = it_left->second;
      auto const& polynomials_right = it_right->second;
      periodic.emplace_hint(
          periodic.cend(),
          ωl,
          typename Result::Polynomials{
              /*sin=*/polynomials_left.sin - polynomials_right.sin,
              /*cos=*/polynomials_left.cos - polynomials_right.cos});
      ++it_left;
      ++it_right;
    }
  }
  return {aperiodic, std::move(periodic)};
}

template<typename Scalar, typename Value, int degree_,
         template<typename, typename, int> class Evaluator>
PoissonSeries<Product<Scalar, Value>, degree_, Evaluator>
operator*(Scalar const& left,
          PoissonSeries<Value, degree_, Evaluator> const& right) {
  using Result = PoissonSeries<Product<Scalar, Value>, degree_, Evaluator>;
  typename Result::PolynomialsByAngularFrequency periodic;
  auto const aperiodic = left * right.aperiodic_;
  for (auto const& [ω, polynomials] : right.periodic_) {
    periodic.emplace_hint(
        periodic.cend(),
        ω,
        typename Result::Polynomials{/*sin=*/left * polynomials.sin,
                                     /*cos=*/left * polynomials.cos});
  }
  return {aperiodic, std::move(periodic)};
}

template<typename Scalar, typename Value, int degree_,
         template<typename, typename, int> class Evaluator>
PoissonSeries<Product<Value, Scalar>, degree_, Evaluator>
operator*(PoissonSeries<Value, degree_, Evaluator> const& left,
          Scalar const& right) {
  using Result = PoissonSeries<Product<Scalar, Value>, degree_, Evaluator>;
  typename Result::PolynomialsByAngularFrequency periodic;
  auto const aperiodic = left.aperiodic_ * right;
  for (auto const& [ω, polynomials] : left.periodic_) {
    periodic.emplace_hint(
        periodic.cend(),
        ω,
        typename Result::Polynomials{/*sin=*/polynomials.sin * right,
                                     /*cos=*/polynomials.cos * right});
  }
  return {aperiodic, std::move(periodic)};
}

template<typename Scalar, typename Value, int degree_,
         template<typename, typename, int> class Evaluator>
PoissonSeries<Quotient<Value, Scalar>, degree_, Evaluator>
operator/(PoissonSeries<Value, degree_, Evaluator> const& left,
          Scalar const& right) {
  using Result = PoissonSeries<Product<Scalar, Value>, degree_, Evaluator>;
  typename Result::PolynomialsByAngularFrequency periodic;
  auto const aperiodic = left.aperiodic_ / right;
  for (auto const& [ω, polynomials] : left.periodic_) {
    periodic.emplace_hint(
        periodic.cend(),
        ω,
        typename Result::Polynomials{/*sin=*/polynomials.sin / right,
                                     /*cos=*/polynomials.cos / right});
  }
  return {aperiodic, std::move(periodic)};
}

template<typename LValue, typename RValue,
         int ldegree_, int rdegree_,
         template<typename, typename, int> class Evaluator>
PoissonSeries<Product<LValue, RValue>, ldegree_ + rdegree_, Evaluator>
operator*(PoissonSeries<LValue, ldegree_, Evaluator> const& left,
          PoissonSeries<RValue, rdegree_, Evaluator> const& right) {
  using Result =
      PoissonSeries<Product<LValue, RValue>, ldegree_ + rdegree_, Evaluator>;

  // Compute all the individual terms using elementary trigonometric identities
  // and put them in a multimap, because the same frequency may appear multiple
  // times.
  std::multimap<AngularFrequency, typename Result::Polynomials> terms;
  auto const aperiodic = left.aperiodic_ * right.aperiodic_;
  for (auto const& [ω, polynomials] : left.periodic_) {
    terms.emplace(ω,
                  typename Result::Polynomials{
                      /*sin=*/polynomials.sin * right.aperiodic_,
                      /*cos=*/polynomials.cos * right.aperiodic_});
  }
  for (auto const& [ω, polynomials] : right.periodic_) {
    terms.emplace(ω,
                  typename Result::Polynomials{
                      /*sin=*/left.aperiodic_ * polynomials.sin,
                      /*cos=*/left.aperiodic_ * polynomials.cos});
  }
  for (auto const& [ωl, polynomials_left] : left.periodic_) {
    for (auto const& [ωr, polynomials_right] : right.periodic_) {
      terms.emplace(ωl - ωr,
                    typename Result::Polynomials{
                        /*sin=*/(-polynomials_left.cos * polynomials_right.sin +
                                 polynomials_left.sin * polynomials_right.cos) /
                            2,
                        /*cos=*/(polynomials_left.sin * polynomials_right.sin +
                                 polynomials_left.cos * polynomials_right.cos) /
                            2});
      terms.emplace(ωl + ωr,
                    typename Result::Polynomials{
                        /*sin=*/(polynomials_left.cos * polynomials_right.sin +
                                 polynomials_left.sin * polynomials_right.cos) /
                            2,
                        /*cos=*/(-polynomials_left.sin * polynomials_right.sin +
                                 polynomials_left.cos * polynomials_right.cos) /
                            2});
    }
  }

  // Now group the terms together by frequency.
  typename Result::PolynomialsByAngularFrequency periodic;
  std::optional<AngularFrequency> previous_ω;
  for (auto it = terms.cbegin(); it != terms.cend(); ++it) {
    auto const& ω = it->first;
    auto const& polynomials = it->second;
    if (previous_ω.has_value() && previous_ω.value() == ω) {
      auto& previous_polynomials = periodic.rbegin()->second;
      previous_polynomials.sin += polynomials.sin;
      previous_polynomials.cos += polynomials.cos;
    } else {
      periodic.insert(*it);
    }
    previous_ω = ω;
  }

  return {aperiodic, std::move(periodic)};
}

template<typename LValue, typename RValue,
         int ldegree_, int rdegree_, int wdegree_,
         template<typename, typename, int> class Evaluator>
Primitive<Product<LValue, RValue>, Time>
Dot(PoissonSeries<LValue, ldegree_, Evaluator> const& left,
    PoissonSeries<RValue, rdegree_, Evaluator> const& right,
    PoissonSeries<double, wdegree_, Evaluator> const& weight,
    Instant const& t_min,
    Instant const& t_max) {
  auto const integrand = left * right * weight;
  auto const primitive = integrand.Primitive();
  return primitive.Evaluate(t_max) - primitive.Evaluate(t_min);
}

template<typename Value, int degree_,
         template<typename, typename, int> class Evaluator>
PiecewisePoissonSeries<Value, degree_, Evaluator>::PiecewisePoissonSeries(
    Interval<Instant> const& interval,
    Series const& series)
    : bounds_({interval.min, interval.max}),
      series_(/*count=*/1, series) {
  CHECK_LT(Time{}, interval.measure());
}

template<typename Value, int degree_,
         template<typename, typename, int> class Evaluator>
void PiecewisePoissonSeries<Value, degree_, Evaluator>::Append(
    Interval<Instant> const& interval,
    Series const& series) {
  CHECK_LT(Time{}, interval.measure());
  CHECK_EQ(bounds_.back(), interval.min);
  bounds_.push_back(interval.max);
  series_.push_back(series);
}

template<typename Value, int degree_,
         template<typename, typename, int> class Evaluator>
Instant PiecewisePoissonSeries<Value, degree_, Evaluator>::t_min() const {
  return bounds_.front();
}

template<typename Value, int degree_,
         template<typename, typename, int> class Evaluator>
Instant PiecewisePoissonSeries<Value, degree_, Evaluator>::t_max() const {
  return bounds_.back();
}

template<typename Value, int degree_,
         template<typename, typename, int> class Evaluator>
Value PiecewisePoissonSeries<Value, degree_, Evaluator>::Evaluate(
    Instant const& t) const {
  // If t is an element of bounds_, the returned iterator points to the next
  // element.  Otherwise it points to the upper bound of the interval to which
  // t belongs.
  auto const it = std::upper_bound(bounds_.cbegin(), bounds_.cend(), t);
  CHECK(it != bounds_.cbegin())
      << "Unexpected result looking up " << t << " in "
      << bounds_.front() << " .. " << bounds_.back();
  CHECK(it != bounds_.cend())
      << t << " is outside of " << bounds_.front() << " .. " << bounds_.back();
  return series_[it - bounds_.cbegin() - 1].Evaluate(t);
}

template<typename Value, int degree_,
         template<typename, typename, int> class Evaluator>
PiecewisePoissonSeries<Value, degree_, Evaluator>::PiecewisePoissonSeries(
    std::vector<Instant> const& bounds,
    std::vector<PoissonSeries<Value, degree_, Evaluator>> const& series)
    : bounds_(bounds),
      series_(series) {}

template<typename Value, int rdegree_,
         template<typename, typename, int> class Evaluator>
PiecewisePoissonSeries<Value, rdegree_, Evaluator> operator+(
    PiecewisePoissonSeries<Value, rdegree_, Evaluator> const& right) {
  return right;
}

template<typename Value, int rdegree_,
         template<typename, typename, int> class Evaluator>
PiecewisePoissonSeries<Value, rdegree_, Evaluator> operator-(
    PiecewisePoissonSeries<Value, rdegree_, Evaluator> const& right) {
  using Result = PiecewisePoissonSeries<Value, rdegree_, Evaluator>;
  std::vector<typename Result::Series> series;
  for (int i = 0; i < right.series_.size(); ++i) {
    series.push_back(-right.series_[i]);
  }
  return Result(right.bounds_, series);
}

template<typename Scalar, typename Value, int degree_,
         template<typename, typename, int> class Evaluator>
PiecewisePoissonSeries<Product<Scalar, Value>, degree_, Evaluator> operator*(
    Scalar const& left,
    PiecewisePoissonSeries<Value, degree_, Evaluator> const& right) {
  using Result =
      PiecewisePoissonSeries<Product<Scalar, Value>, degree_, Evaluator>;
  std::vector<typename Result::Series> series;
  for (int i = 0; i < right.series_.size(); ++i) {
    series.push_back(left * right.series_[i]);
  }
  return Result(right.bounds_, series);
}

template<typename Scalar, typename Value, int degree_,
         template<typename, typename, int> class Evaluator>
PiecewisePoissonSeries<Product<Value, Scalar>, degree_, Evaluator> operator*(
    PiecewisePoissonSeries<Value, degree_, Evaluator> const& left,
    Scalar const& right) {
  using Result =
      PiecewisePoissonSeries<Product<Value, Scalar>, degree_, Evaluator>;
  std::vector<typename Result::Series> series;
  for (int i = 0; i < left.series_.size(); ++i) {
    series.push_back(left.series_[i] * right);
  }
  return Result(left.bounds_, series);
}

template<typename Scalar, typename Value, int degree_,
         template<typename, typename, int> class Evaluator>
PiecewisePoissonSeries<Quotient<Value, Scalar>, degree_, Evaluator> operator/(
    PiecewisePoissonSeries<Value, degree_, Evaluator> const& left,
    Scalar const& right) {
  using Result =
      PiecewisePoissonSeries<Quotient<Value, Scalar>, degree_, Evaluator>;
  std::vector<typename Result::Series> series;
  for (int i = 0; i < left.series_.size(); ++i) {
    series.push_back(left.series_[i] / right);
  }
  return Result(left.bounds_, series);
}

template<typename Value, int ldegree_, int rdegree_,
         template<typename, typename, int> class Evaluator>
PiecewisePoissonSeries<Value, std::max(ldegree_, rdegree_), Evaluator>
operator+(PoissonSeries<Value, ldegree_, Evaluator> const& left,
          PiecewisePoissonSeries<Value, rdegree_, Evaluator> const& right) {
  using Result =
      PiecewisePoissonSeries<Value, std::max(ldegree_, rdegree_), Evaluator>;
  std::vector<typename Result::Series> series;
  for (int i = 0; i < right.series_.size(); ++i) {
    series.push_back(left + right.series_[i]);
  }
  return Result(right.bounds_, series);
}

template<typename Value, int ldegree_, int rdegree_,
         template<typename, typename, int> class Evaluator>
PiecewisePoissonSeries<Value, std::max(ldegree_, rdegree_), Evaluator>
operator+(PiecewisePoissonSeries<Value, ldegree_, Evaluator> const& left,
          PoissonSeries<Value, rdegree_, Evaluator> const& right) {
  using Result =
      PiecewisePoissonSeries<Value, std::max(ldegree_, rdegree_), Evaluator>;
  std::vector<typename Result::Series> series;
  for (int i = 0; i < left.series_.size(); ++i) {
    series.push_back(left.series_[i] + right);
  }
  return Result(left.bounds_, series);
}

template<typename Value, int ldegree_, int rdegree_,
         template<typename, typename, int> class Evaluator>
PiecewisePoissonSeries<Value, std::max(ldegree_, rdegree_), Evaluator>
operator-(PoissonSeries<Value, ldegree_, Evaluator> const& left,
          PiecewisePoissonSeries<Value, rdegree_, Evaluator> const& right) {
  using Result =
      PiecewisePoissonSeries<Value, std::max(ldegree_, rdegree_), Evaluator>;
  std::vector<typename Result::Series> series;
  for (int i = 0; i < right.series_.size(); ++i) {
    series.push_back(left - right.series_[i]);
  }
  return Result(right.bounds_, series);
}

template<typename Value, int ldegree_, int rdegree_,
         template<typename, typename, int> class Evaluator>
PiecewisePoissonSeries<Value, std::max(ldegree_, rdegree_), Evaluator>
operator-(PiecewisePoissonSeries<Value, ldegree_, Evaluator> const& left,
          PoissonSeries<Value, rdegree_, Evaluator> const& right) {
  using Result =
      PiecewisePoissonSeries<Value, std::max(ldegree_, rdegree_), Evaluator>;
  std::vector<typename Result::Series> series;
  for (int i = 0; i < left.series_.size(); ++i) {
    series.push_back(left.series_[i] - right);
  }
  return Result(left.bounds_, series);
}

template<typename LValue, typename RValue, int ldegree_, int rdegree_,
         template<typename, typename, int> class Evaluator>
PiecewisePoissonSeries<Product<LValue, RValue>, ldegree_ + rdegree_, Evaluator>
operator*(PoissonSeries<LValue, ldegree_, Evaluator> const& left,
          PiecewisePoissonSeries<RValue, rdegree_, Evaluator> const& right) {
  using Result = PiecewisePoissonSeries<Product<LValue, RValue>,
                                        ldegree_ + rdegree_,
                                        Evaluator>;
  std::vector<typename Result::Series> series;
  for (int i = 0; i < right.series_.size(); ++i) {
    series.push_back(left * right.series_[i]);
  }
  return Result(right.bounds_, series);
}

template<typename LValue, typename RValue, int ldegree_, int rdegree_,
         template<typename, typename, int> class Evaluator>
PiecewisePoissonSeries<Product<LValue, RValue>, ldegree_ + rdegree_, Evaluator>
operator*(PiecewisePoissonSeries<LValue, ldegree_, Evaluator> const& left,
          PoissonSeries<RValue, rdegree_, Evaluator> const& right) {
  using Result = PiecewisePoissonSeries<Product<LValue, RValue>,
                                        ldegree_ + rdegree_,
                                        Evaluator>;
  std::vector<typename Result::Series> series;
  for (int i = 0; i < left.series_.size(); ++i) {
    series.push_back(left.series_[i] * right);
  }
  return Result(left.bounds_, series);
}

template<typename LValue, typename RValue,
         int ldegree_, int rdegree_, int wdegree_,
         template<typename, typename, int> class Evaluator>
Primitive<Product<LValue, RValue>, Time> Dot(
    PoissonSeries<LValue, ldegree_, Evaluator> const& left,
    PiecewisePoissonSeries<RValue, rdegree_, Evaluator> const& right,
    PoissonSeries<double, wdegree_, Evaluator> const& weight) {
  using Result = Primitive<Product<LValue, RValue>, Time>;
  Result result;
  std::vector<typename Result::Series> series;
  for (int i = 0; i < right.series_.size(); ++i) {
    //TODO(phl):correct use of weight?
    result += Dot(left,
                  right.series_[i],
                  weight,
                  right.bounds_[i],
                  right.bounds_[i + 1]);
  }
  return result;
}

template<typename LValue, typename RValue,
         int ldegree_, int rdegree_, int wdegree_,
         template<typename, typename, int> class Evaluator>
Primitive<Product<LValue, RValue>, Time> Dot(
    PiecewisePoissonSeries<LValue, ldegree_, Evaluator> const& left,
    PoissonSeries<RValue, rdegree_, Evaluator> const& right,
    PoissonSeries<double, wdegree_, Evaluator> const& weight) {
  using Result = Primitive<Product<LValue, RValue>, Time>;
  Result result;
  std::vector<typename Result::Series> series;
  for (int i = 0; i < left.series_.size(); ++i) {
    //TODO(phl):correct use of weight?
    result += Dot(left.series_[i],
                  right,
                  weight,
                  left.bounds_[i],
                  left.bounds_[i + 1]);
  }
  return result;
}

}  // namespace internal_poisson_series
}  // namespace numerics
}  // namespace principia
