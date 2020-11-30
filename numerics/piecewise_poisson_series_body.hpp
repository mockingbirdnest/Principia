
#pragma once

#include "numerics/piecewise_poisson_series.hpp"

#include <algorithm>
#include <vector>

#include "quantities/elementary_functions.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/numerics.hpp"

namespace principia {
namespace numerics {
namespace internal_piecewise_poisson_series {

using quantities::Cos;
using quantities::Sin;
using quantities::si::Radian;

#define USE_GAUSS 0
//TODO(phl):Cleanup
// The minimum value of the max_point parameter passed to Clenshaw-Curtis
// integration, irrespective of the frequencies of the argument function.
constexpr int clenshaw_curtis_min_points_overall = 65;

// The maximum number of points use in Clenshaw-Curtis integration for each
// period of the highest frequency of the argument function.
constexpr int clenshaw_curtis_point_per_period = 5;

// The desired relative error on Clenshaw-Curtis integration, as determined by
// two successive computations with increasing number of points.
constexpr double clenshaw_curtis_relative_error = 0x1p-32;

template<typename Value,
         int aperiodic_degree_, int periodic_degree_,
         template<typename, typename, int> class Evaluator>
PiecewisePoissonSeries<Value, aperiodic_degree_, periodic_degree_, Evaluator>::
PiecewisePoissonSeries(Interval<Instant> const& interval,
                        Series const& series)
    : bounds_({interval.min, interval.max}),
      series_(/*count=*/1, series) {
  CHECK_LT(Time{}, interval.measure());
}

template<typename Value,
         int aperiodic_degree_, int periodic_degree_,
         template<typename, typename, int> class Evaluator>
void
PiecewisePoissonSeries<Value, aperiodic_degree_, periodic_degree_, Evaluator>::
Append(Interval<Instant> const& interval,
       Series const& series) {
  CHECK_LT(Time{}, interval.measure());
  CHECK_EQ(bounds_.back(), interval.min);
  bounds_.push_back(interval.max);
  series_.push_back(series);
}

template<typename Value,
         int aperiodic_degree_, int periodic_degree_,
         template<typename, typename, int> class Evaluator>
Instant
PiecewisePoissonSeries<Value, aperiodic_degree_, periodic_degree_, Evaluator>::
t_min() const {
  return bounds_.front();
}

template<typename Value,
         int aperiodic_degree_, int periodic_degree_,
         template<typename, typename, int> class Evaluator>
Instant
PiecewisePoissonSeries<Value, aperiodic_degree_, periodic_degree_, Evaluator>::
t_max() const {
  return bounds_.back();
}

template<typename Value,
         int aperiodic_degree_, int periodic_degree_,
         template<typename, typename, int> class Evaluator>
Value
PiecewisePoissonSeries<Value, aperiodic_degree_, periodic_degree_, Evaluator>::
operator()(Instant const& t) const {
  Value const addend = EvaluateAddend(t);
  if (t == bounds_.back()) {
    return series_.back()(t) + addend;
  }

  // If t is an element of bounds_, the returned iterator points to the next
  // element.  Otherwise it points to the upper bound of the interval to which
  // t belongs.
  auto const it = std::upper_bound(bounds_.cbegin(), bounds_.cend(), t);
  CHECK(it != bounds_.cbegin())
      << "Unexpected result looking up " << t << " in "
      << bounds_.front() << " .. " << bounds_.back();
  CHECK(it != bounds_.cend())
      << t << " is outside of " << bounds_.front() << " .. " << bounds_.back();
  return series_[it - bounds_.cbegin() - 1](t) + addend;
}

template<typename Value,
         int aperiodic_degree_, int periodic_degree_,
         template<typename, typename, int> class Evaluator>
auto
PiecewisePoissonSeries<Value, aperiodic_degree_, periodic_degree_, Evaluator>::
FourierTransform() const -> Spectrum {
  // TODO(egg): consider pre-evaluating |*this| at all points used by the
  // Gaussian quadratures, removing the lifetime requirement on |*this| and
  // potentially speeding up repeated evaluations of the Fourier transform.
  return [this](AngularFrequency const& ω) {
    Interval<Instant> const time_domain{t_min(), t_max()};
    Instant const t0 = time_domain.midpoint();
    Primitive<Complexification<Value>, Instant> integral;
    for (int k = 0; k < series_.size(); ++k) {
      integral +=
          quadrature::GaussLegendre<std::max(1, (aperiodic_degree_ + 1) / 2)>(
              [this, &f = series_[k], t0, ω](
                  Instant const& t) -> Complexification<Value> {
                return (f(t) + EvaluateAddend(t)) *
                       Complexification<double>{Cos(ω * (t - t0)),
                                                -Sin(ω * (t - t0))};
              },
              bounds_[k],
              bounds_[k + 1]);
    }
    return integral;
  };
}

template<typename Value,
         int aperiodic_degree_, int periodic_degree_,
         template<typename, typename, int> class Evaluator>
template<int aperiodic_rdegree, int periodic_rdegree>
PiecewisePoissonSeries<Value, aperiodic_degree_, periodic_degree_, Evaluator>&
PiecewisePoissonSeries<Value, aperiodic_degree_, periodic_degree_, Evaluator>::
operator+=(PoissonSeries<Value,
                         aperiodic_rdegree, periodic_rdegree,
                         Evaluator> const& right) {
  static_assert(aperiodic_rdegree <= aperiodic_degree_);
  static_assert(periodic_rdegree <= periodic_degree_);
  if (addend_.has_value()) {
    addend_.value() += right;
  } else {
    addend_ = Series(right);
  }
  return *this;
}

template<typename Value,
         int aperiodic_degree_, int periodic_degree_,
         template<typename, typename, int> class Evaluator>
template<int aperiodic_rdegree, int periodic_rdegree>
PiecewisePoissonSeries<Value, aperiodic_degree_, periodic_degree_, Evaluator>&
PiecewisePoissonSeries<Value, aperiodic_degree_, periodic_degree_, Evaluator>::
operator-=(PoissonSeries<Value,
                         aperiodic_rdegree, periodic_rdegree,
                         Evaluator> const& right) {
  static_assert(aperiodic_rdegree <= aperiodic_degree_);
  static_assert(periodic_rdegree <= periodic_degree_);
  if (addend_.has_value()) {
    addend_.value() -= right;
  } else {
    addend_ = Series(-right);
  }
  return *this;
}

template<typename Value,
         int aperiodic_degree_, int periodic_degree_,
         template<typename, typename, int> class Evaluator>
void
PiecewisePoissonSeries<Value, aperiodic_degree_, periodic_degree_, Evaluator>::
WriteToMessage(not_null<serialization::PiecewisePoissonSeries*> message) const {
  for (Instant const& bound : bounds_) {
    bound.WriteToMessage(message->add_bounds());
  }
  for (Series const& series : series_) {
    series.WriteToMessage(message->add_series());
  }
  if (addend_.has_value()) {
    addend_->WriteToMessage(message->mutable_addend());
  }
}

template<typename Value,
         int aperiodic_degree_, int periodic_degree_,
         template<typename, typename, int> class Evaluator>
PiecewisePoissonSeries<Value, aperiodic_degree_, periodic_degree_, Evaluator>
PiecewisePoissonSeries<Value, aperiodic_degree_, periodic_degree_, Evaluator>::
ReadFromMessage(serialization::PiecewisePoissonSeries const& message) {
  CHECK_NE(0, message.series_size());
  CHECK_EQ(message.bounds_size(), message.series_size() + 1);
  Interval<Instant> const first_interval{
      Instant::ReadFromMessage(message.bounds(0)),
      Instant::ReadFromMessage(message.bounds(1))};
  PiecewisePoissonSeries series(first_interval,
                                Series::ReadFromMessage(message.series(0)));
  for (int i = 1; i < message.series_size(); ++i) {
    Interval<Instant> const interval{
        Instant::ReadFromMessage(message.bounds(i)),
        Instant::ReadFromMessage(message.bounds(i + 1))};
    series.Append(interval, Series::ReadFromMessage(message.series(i)));
  }
  if (message.has_addend()) {
    series.addend_ = Series::ReadFromMessage(message.addend());
  }
  return series;
}

template<typename Value,
         int aperiodic_degree_, int periodic_degree_,
         template<typename, typename, int> class Evaluator>
PiecewisePoissonSeries<Value, aperiodic_degree_, periodic_degree_, Evaluator>::
PiecewisePoissonSeries(std::vector<Instant> const& bounds,
                       std::vector<PoissonSeries<Value,
                                                 aperiodic_degree_,
                                                 periodic_degree_,
                                                 Evaluator>> const& series,
                       std::optional<Series> const& addend)
    : bounds_(bounds),
      series_(series),
      addend_(addend) {}

template<typename Value,
         int aperiodic_degree_, int periodic_degree_,
         template<typename, typename, int> class Evaluator>
Value
PiecewisePoissonSeries<Value, aperiodic_degree_, periodic_degree_, Evaluator>::
EvaluateAddend(Instant const& t) const {
  return addend_.has_value() ? addend_.value()(t) : Value{};
}

template<typename Value,
         int aperiodic_rdegree, int periodic_rdegree,
         template<typename, typename, int> class Evaluator>
PiecewisePoissonSeries<Value, aperiodic_rdegree, periodic_rdegree, Evaluator>
operator+(PiecewisePoissonSeries<Value,
                                 aperiodic_rdegree, periodic_rdegree,
                                 Evaluator> const& right) {
  return right;
}

template<typename Value,
         int aperiodic_rdegree, int periodic_rdegree,
         template<typename, typename, int> class Evaluator>
PiecewisePoissonSeries<Value, aperiodic_rdegree, periodic_rdegree, Evaluator>
operator-(PiecewisePoissonSeries<Value,
                                 aperiodic_rdegree, periodic_rdegree,
                                 Evaluator> const& right) {
  using Result = PiecewisePoissonSeries<Value,
                                        aperiodic_rdegree, periodic_rdegree,
                                        Evaluator>;
  std::vector<typename Result::Series> series;
  series.reserve(right.series_.size());
  for (int i = 0; i < right.series_.size(); ++i) {
    series.push_back(-right.series_[i]);
  }
  std::optional<typename Result::Series> addend;
  if (right.addend_.has_value()) {
    addend = -right.addend_.value();
  }
  return Result(right.bounds_, series, addend);
}

template<typename Scalar,
         typename Value,
         int aperiodic_rdegree, int periodic_rdegree,
         template<typename, typename, int> class Evaluator>
PiecewisePoissonSeries<Product<Scalar, Value>,
                       aperiodic_rdegree, periodic_rdegree,
                       Evaluator>
operator*(Scalar const& left,
          PiecewisePoissonSeries<Value,
                                 aperiodic_rdegree, periodic_rdegree,
                                 Evaluator> const& right) {
  using Result = PiecewisePoissonSeries<Product<Scalar, Value>,
                                        aperiodic_rdegree, periodic_rdegree,
                                        Evaluator>;
  std::vector<typename Result::Series> series;
  series.reserve(right.series_.size());
  for (int i = 0; i < right.series_.size(); ++i) {
    series.push_back(left * right.series_[i]);
  }
  std::optional<typename Result::Series> addend;
  if (right.addend_.has_value()) {
    addend = left * right.addend_.value();
  }
  return Result(right.bounds_, series, addend);
}

template<typename Scalar,
         typename Value,
         int aperiodic_ldegree, int periodic_ldegree,
         template<typename, typename, int> class Evaluator>
PiecewisePoissonSeries<Product<Value, Scalar>,
                       aperiodic_ldegree, periodic_ldegree,
                       Evaluator>
operator*(PiecewisePoissonSeries<Value,
                                 aperiodic_ldegree, periodic_ldegree,
                                 Evaluator> const& left,
          Scalar const& right) {
  using Result = PiecewisePoissonSeries<Product<Value, Scalar>,
                                        aperiodic_ldegree, periodic_ldegree,
                                        Evaluator>;
  std::vector<typename Result::Series> series;
  series.reserve(left.series_.size());
  for (int i = 0; i < left.series_.size(); ++i) {
    series.push_back(left.series_[i] * right);
  }
  std::optional<typename Result::Series> addend;
  if (left.addend_.has_value()) {
    addend = left.addend_.value() * right;
  }
  return Result(left.bounds_, series, addend);
}

template<typename Scalar,
         typename Value,
         int aperiodic_ldegree, int periodic_ldegree,
         template<typename, typename, int> class Evaluator>
PiecewisePoissonSeries<Quotient<Value, Scalar>,
                       aperiodic_ldegree, periodic_ldegree,
                       Evaluator>
operator/(PiecewisePoissonSeries<Value,
                                 aperiodic_ldegree, periodic_ldegree,
                                 Evaluator> const& left,
          Scalar const& right) {
  using Result = PiecewisePoissonSeries<Quotient<Value, Scalar>,
                                        aperiodic_ldegree, periodic_ldegree,
                                        Evaluator>;
  std::vector<typename Result::Series> series;
  series.reserve(left.series_.size());
  for (int i = 0; i < left.series_.size(); ++i) {
    series.push_back(left.series_[i] / right);
  }
  std::optional<typename Result::Series> addend;
  if (left.addend_.has_value()) {
    addend = left.addend_.value() / right;
  }
  return Result(left.bounds_, series, addend);
}

template<typename Value,
         int aperiodic_ldegree, int periodic_ldegree,
         int aperiodic_rdegree, int periodic_rdegree,
         template<typename, typename, int> class Evaluator>
PiecewisePoissonSeries<Value,
                       PRINCIPIA_MAX(aperiodic_ldegree, aperiodic_rdegree),
                       PRINCIPIA_MAX(periodic_ldegree, periodic_rdegree),
                       Evaluator>
operator+(PoissonSeries<Value,
                        aperiodic_ldegree, periodic_ldegree,
                        Evaluator> const& left,
          PiecewisePoissonSeries<Value,
                                 aperiodic_rdegree, periodic_rdegree,
                                 Evaluator> const& right) {
  using Result =
      PiecewisePoissonSeries<Value,
                             std::max(aperiodic_ldegree, aperiodic_rdegree),
                             std::max(periodic_ldegree, periodic_rdegree),
                             Evaluator>;
  std::optional<typename Result::Series> addend;
  if (right.addend_.has_value()) {
    addend = left + right.addend_.value();
  } else {
    addend = typename Result::Series(left);
  }
  return Result(right.bounds_, right.series_, addend);
}

template<typename Value,
         int aperiodic_ldegree, int periodic_ldegree,
         int aperiodic_rdegree, int periodic_rdegree,
         template<typename, typename, int> class Evaluator>
PiecewisePoissonSeries<Value,
                       PRINCIPIA_MAX(aperiodic_ldegree, aperiodic_rdegree),
                       PRINCIPIA_MAX(periodic_ldegree, periodic_rdegree),
                       Evaluator>
operator+(PiecewisePoissonSeries<Value,
                                 aperiodic_ldegree, periodic_ldegree,
                                 Evaluator> const& left,
          PoissonSeries<Value,
                        aperiodic_rdegree, periodic_rdegree,
                        Evaluator> const& right) {
  using Result =
      PiecewisePoissonSeries<Value,
                             std::max(aperiodic_ldegree, aperiodic_rdegree),
                             std::max(periodic_ldegree, periodic_rdegree),
                             Evaluator>;
  std::optional<typename Result::Series> addend;
  if (left.addend_.has_value()) {
    addend = left.addend_.value() + right;
  } else {
    addend = typename Result::Series(right);
  }
  return Result(left.bounds_, left.series_, addend);
}

template<typename Value,
         int aperiodic_ldegree, int periodic_ldegree,
         int aperiodic_rdegree, int periodic_rdegree,
         template<typename, typename, int> class Evaluator>
PiecewisePoissonSeries<Value,
                       PRINCIPIA_MAX(aperiodic_ldegree, aperiodic_rdegree),
                       PRINCIPIA_MAX(periodic_ldegree, periodic_rdegree),
                       Evaluator>
operator-(PoissonSeries<Value,
                        aperiodic_ldegree, periodic_ldegree,
                        Evaluator> const& left,
          PiecewisePoissonSeries<Value,
                                 aperiodic_rdegree, periodic_rdegree,
                                 Evaluator> const& right) {
  using Result =
      PiecewisePoissonSeries<Value,
                             std::max(aperiodic_ldegree, aperiodic_rdegree),
                             std::max(periodic_ldegree, periodic_rdegree),
                             Evaluator>;
  std::vector<typename Result::Series> series;
  series.reserve(right.series_.size());
  for (int i = 0; i < right.series_.size(); ++i) {
    series.push_back(-right.series_[i]);
  }
  std::optional<typename Result::Series> addend;
  if (right.addend_.has_value()) {
    addend = left - right.addend_.value();
  } else {
    addend = typename Result::Series(left);
  }
  return Result(right.bounds_, series, addend);
}

template<typename Value,
         int aperiodic_ldegree, int periodic_ldegree,
         int aperiodic_rdegree, int periodic_rdegree,
         template<typename, typename, int> class Evaluator>
PiecewisePoissonSeries<Value,
                       PRINCIPIA_MAX(aperiodic_ldegree, aperiodic_rdegree),
                       PRINCIPIA_MAX(periodic_ldegree, periodic_rdegree),
                       Evaluator>
operator-(PiecewisePoissonSeries<Value,
                                 aperiodic_ldegree, periodic_ldegree,
                                 Evaluator> const& left,
          PoissonSeries<Value,
                        aperiodic_rdegree, periodic_rdegree,
                        Evaluator> const& right) {
  using Result =
      PiecewisePoissonSeries<Value,
                             std::max(aperiodic_ldegree, aperiodic_rdegree),
                             std::max(periodic_ldegree, periodic_rdegree),
                             Evaluator>;
  std::vector<typename Result::Series> series;
  series.reserve(left.series_.size());
  for (int i = 0; i < left.series_.size(); ++i) {
    series.push_back(typename Result::Series(left.series_[i]));
  }
  std::optional<typename Result::Series> addend;
  if (left.addend_.has_value()) {
    addend = left.addend_.value() - right;
  } else {
    addend = typename Result::Series(-right);
  }
  return Result(left.bounds_, series, addend);
}

template<typename LValue, typename RValue,
         int aperiodic_ldegree, int periodic_ldegree,
         int aperiodic_rdegree, int periodic_rdegree,
         template<typename, typename, int> class Evaluator>
PiecewisePoissonSeries<Product<LValue, RValue>,
                       PRINCIPIA_MAX4(aperiodic_ldegree + aperiodic_rdegree,
                                      aperiodic_ldegree + periodic_rdegree,
                                      periodic_ldegree + aperiodic_rdegree,
                                      periodic_ldegree + periodic_rdegree),
                       PRINCIPIA_MAX3(aperiodic_ldegree + periodic_rdegree,
                                      periodic_ldegree + aperiodic_rdegree,
                                      periodic_ldegree + periodic_rdegree),
                       Evaluator>
operator*(PoissonSeries<LValue,
                        aperiodic_ldegree, periodic_ldegree,
                        Evaluator> const& left,
          PiecewisePoissonSeries<RValue,
                                 aperiodic_rdegree, periodic_rdegree,
                                 Evaluator> const& right) {
  using Result =
      PiecewisePoissonSeries<Product<LValue, RValue>,
                             std::max({aperiodic_ldegree + aperiodic_rdegree,
                                       aperiodic_ldegree + periodic_rdegree,
                                       periodic_ldegree + aperiodic_rdegree,
                                       periodic_ldegree + periodic_rdegree}),
                             std::max({aperiodic_ldegree + periodic_rdegree,
                                       periodic_ldegree + aperiodic_rdegree,
                                       periodic_ldegree + periodic_rdegree}),
                             Evaluator>;
  std::vector<typename Result::Series> series;
  series.reserve(right.series_.size());
  for (int i = 0; i < right.series_.size(); ++i) {
    Instant const origin = right.series_[i].origin();
    series.push_back(left.AtOrigin(origin) * right.series_[i]);
  }
  std::optional<typename Result::Series> addend;
  if (right.addend_.has_value()) {
    addend = left * right.addend_.value();
  }
  return Result(right.bounds_, series, addend);
}

template<typename LValue, typename RValue,
         int aperiodic_ldegree, int periodic_ldegree,
         int aperiodic_rdegree, int periodic_rdegree,
         template<typename, typename, int> class Evaluator>
PiecewisePoissonSeries<Product<LValue, RValue>,
                       PRINCIPIA_MAX4(aperiodic_ldegree + aperiodic_rdegree,
                                      aperiodic_ldegree + periodic_rdegree,
                                      periodic_ldegree + aperiodic_rdegree,
                                      periodic_ldegree + periodic_rdegree),
                       PRINCIPIA_MAX3(aperiodic_ldegree + periodic_rdegree,
                                      periodic_ldegree + aperiodic_rdegree,
                                      periodic_ldegree + periodic_rdegree),
                       Evaluator>
operator*(PiecewisePoissonSeries<LValue,
                                 aperiodic_ldegree, periodic_ldegree,
                                 Evaluator> const& left,
          PoissonSeries<RValue,
                        aperiodic_rdegree, periodic_rdegree,
                        Evaluator> const& right) {
  using Result =
      PiecewisePoissonSeries<Product<LValue, RValue>,
                             std::max({aperiodic_ldegree + aperiodic_rdegree,
                                       aperiodic_ldegree + periodic_rdegree,
                                       periodic_ldegree + aperiodic_rdegree,
                                       periodic_ldegree + periodic_rdegree}),
                             std::max({aperiodic_ldegree + periodic_rdegree,
                                       periodic_ldegree + aperiodic_rdegree,
                                       periodic_ldegree + periodic_rdegree}),
                             Evaluator>;
  std::vector<typename Result::Series> series;
  series.reserve(left.series_.size());
  for (int i = 0; i < left.series_.size(); ++i) {
    Instant const origin = left.series_[i].origin();
    series.push_back(left.series_[i] * right.AtOrigin(origin));
  }
  std::optional<typename Result::Series> addend;
  if (left.addend_.has_value()) {
    addend = left.addend_.value() * right;
  }
  return Result(left.bounds_, series, addend);
}

template<typename LValue, typename RValue,
         int aperiodic_ldegree, int periodic_ldegree,
         int aperiodic_rdegree, int periodic_rdegree,
         int aperiodic_wdegree, int periodic_wdegree,
         template<typename, typename, int> class Evaluator,
         int points>
typename Hilbert<LValue, RValue>::InnerProductType
InnerProduct(PoissonSeries<LValue,
                           aperiodic_ldegree, periodic_ldegree,
                           Evaluator> const& left,
             PiecewisePoissonSeries<RValue,
                                    aperiodic_rdegree, periodic_rdegree,
                                    Evaluator> const& right,
             PoissonSeries<double,
                           aperiodic_wdegree, periodic_wdegree,
                           Evaluator> const& weight) {
  return InnerProduct<LValue, RValue,
                      aperiodic_ldegree, periodic_ldegree,
                      aperiodic_rdegree, periodic_rdegree,
                      aperiodic_wdegree, periodic_wdegree,
                      Evaluator,
                      points>(
      left, right, weight, right.t_min(), right.t_max());
}

template<typename LValue, typename RValue,
         int aperiodic_ldegree, int periodic_ldegree,
         int aperiodic_rdegree, int periodic_rdegree,
         int aperiodic_wdegree, int periodic_wdegree,
         template<typename, typename, int> class Evaluator,
         int points>
typename Hilbert<LValue, RValue>::InnerProductType
InnerProduct(PoissonSeries<LValue,
                           aperiodic_ldegree, periodic_ldegree,
                           Evaluator> const& left,
             PiecewisePoissonSeries<RValue,
                                    aperiodic_rdegree, periodic_rdegree,
                                    Evaluator> const& right,
             PoissonSeries<double,
                           aperiodic_wdegree, periodic_wdegree,
                           Evaluator> const& weight,
             Instant const& t_min,
             Instant const& t_max) {
#if USE_GAUSS
  using Result =
      Primitive<typename Hilbert<LValue, RValue>::InnerProductType, Time>;
  Result result{};
  for (int i = 0; i < right.series_.size(); ++i) {
    auto integrand = [i, &left, &right, &weight](Instant const& t) {
      return Hilbert<LValue, RValue>::InnerProduct(
          left(t) * weight(t),
          right.series_[i](t) + right.EvaluateAddend(t));
    };
    auto const integral = quadrature::GaussLegendre<points>(
        integrand, right.bounds_[i], right.bounds_[i + 1]);
    result += integral;
  }
  return result / (t_max - t_min);
#else
  AngularFrequency max_ω =
      left.addend_.has_value() ? right.addend_->max_ω() : AngularFrequency{};
  for (int i = 0; i < right.series_.size(); ++i) {
    max_ω = std::max(max_ω, right.series_[i].max_ω());
  }
  max_ω += left.max_ω() + weight.max_ω();

  std::optional<int> max_points =
      max_ω == AngularFrequency()
          ? std::optional<int>{}
          : std::max(
                clenshaw_curtis_min_points_overall,
                static_cast<int>(clenshaw_curtis_point_per_period *
                                 (t_max - t_min) * max_ω / (2 * π * Radian)));

  auto integrand = [&left, &right, &weight](Instant const& t) {
    return Hilbert<LValue, RValue>::InnerProduct(left(t), right(t)) * weight(t);
  };
  return quadrature::AutomaticClenshawCurtis(
             integrand,
             t_min,
             t_max,
             /*max_relative_error=*/clenshaw_curtis_relative_error,
             max_points) /
         (t_max - t_min);
#endif
}

template<typename LValue, typename RValue,
         int aperiodic_ldegree, int periodic_ldegree,
         int aperiodic_rdegree, int periodic_rdegree,
         int aperiodic_wdegree, int periodic_wdegree,
         template<typename, typename, int> class Evaluator,
         int points>
typename Hilbert<LValue, RValue>::InnerProductType
InnerProduct(PiecewisePoissonSeries<LValue,
                                    aperiodic_ldegree, periodic_ldegree,
                                    Evaluator> const& left,
             PoissonSeries<RValue,
                           aperiodic_rdegree, periodic_rdegree,
                           Evaluator> const& right,
             PoissonSeries<double,
                           aperiodic_wdegree, periodic_wdegree,
                           Evaluator> const& weight) {
  return InnerProduct<LValue, RValue,
                      aperiodic_ldegree, periodic_ldegree,
                      aperiodic_rdegree, periodic_rdegree,
                      aperiodic_wdegree, periodic_wdegree,
                      Evaluator,
                      points>(
       left, right, weight, left.t_min(), left.t_max());
}

template<typename LValue, typename RValue,
         int aperiodic_ldegree, int periodic_ldegree,
         int aperiodic_rdegree, int periodic_rdegree,
         int aperiodic_wdegree, int periodic_wdegree,
         template<typename, typename, int> class Evaluator,
         int points>
typename Hilbert<LValue, RValue>::InnerProductType
InnerProduct(PiecewisePoissonSeries<LValue,
                                    aperiodic_ldegree, periodic_ldegree,
                                    Evaluator> const& left,
             PoissonSeries<RValue,
                           aperiodic_rdegree, periodic_rdegree,
                           Evaluator> const& right,
             PoissonSeries<double,
                           aperiodic_wdegree, periodic_wdegree,
                           Evaluator> const& weight,
             Instant const& t_min,
             Instant const& t_max) {
#if USE_GAUSS
  using Result =
      Primitive<typename Hilbert<LValue, RValue>::InnerProductType, Time>;
  Result result{};
  for (int i = 0; i < left.series_.size(); ++i) {
    auto integrand = [i, &left, &right, &weight](Instant const& t) {
      return Hilbert<LValue, RValue>::InnerProduct(
          left.series_[i](t) + left.EvaluateAddend(t),
          right(t) * weight(t));
    };
    auto const integral = quadrature::GaussLegendre<points>(
        integrand, left.bounds_[i], left.bounds_[i + 1]);
    result += integral;
  }
  if constexpr (points == 15 || points == 10 || points == 7) {
    auto const ip1 = result / (t_max - t_min);
    auto const ip2 =
        InnerProduct<LValue,
                     RValue,
                     aperiodic_ldegree,
                     periodic_ldegree,
                     aperiodic_rdegree,
                     periodic_rdegree,
                     aperiodic_wdegree,
                     periodic_wdegree,
                     Evaluator,
                     2 * points>(left, right, weight, t_min, t_max);
    static double max_relative_error = 0;
    double const relative_error = testing_utilities::RelativeError(ip1, ip2);
    if (relative_error > max_relative_error) {
      max_relative_error = relative_error;
      LOG(ERROR) << max_relative_error;
    }
    return ip2;
  } else {
    return result / (t_max - t_min);
  }
#else
  AngularFrequency max_ω =
      left.addend_.has_value() ? left.addend_->max_ω() : AngularFrequency{};
  for (int i = 0; i < left.series_.size(); ++i) {
    max_ω = std::max(max_ω, left.series_[i].max_ω());
  }
  max_ω += right.max_ω() + weight.max_ω();

  //TODO(phl):Cleanup
  std::optional<int> max_points =
      max_ω == AngularFrequency()
          ? std::optional<int>{}
          : std::max(
                clenshaw_curtis_min_points_overall,
                static_cast<int>(clenshaw_curtis_point_per_period *
                                 (t_max - t_min) * max_ω / (2 * π * Radian)));

  auto integrand = [&left, &right, &weight](Instant const& t) {
    return Hilbert<LValue, RValue>::InnerProduct(left(t), right(t)) * weight(t);
  };
  quadrature::internal_quadrature::do_the_logging = true;
  auto const result = quadrature::AutomaticClenshawCurtis(
             integrand,
             t_min,
             t_max,
             /*max_relative_error=*/clenshaw_curtis_relative_error,
             max_points) /
         (t_max - t_min);
  quadrature::internal_quadrature::do_the_logging = false;
  return result;
#endif
}

}  // namespace internal_piecewise_poisson_series
}  // namespace numerics
}  // namespace principia
