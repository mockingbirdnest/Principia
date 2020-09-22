
#pragma once

#include "numerics/poisson_series.hpp"

#include <algorithm>
#include <functional>
#include <map>
#include <optional>
#include <utility>
#include <vector>

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
using quantities::Primitive;
using quantities::Sin;
using quantities::Time;
using quantities::Variation;
using quantities::si::Metre;
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
  using Result =
      PoissonSeries<typename Hilbert<LValue, RValue>::InnerProductType,
                    ldegree_ + rdegree_,
                    Evaluator>;

  auto aperiodic = product(left.aperiodic_, right.aperiodic_);

  // Compute all the individual terms using elementary trigonometric identities
  // and put them in a vector, because the same frequency may appear multiple
  // times.
  using Term = std::pair<AngularFrequency, typename Result::Polynomials>;
  std::vector<Term> terms;
  terms.reserve(left.periodic_.size() + right.periodic_.size() +
                left.periodic_.size() * right.periodic_.size());
  for (auto const& [ω, polynomials] : left.periodic_) {
    terms.emplace_back(ω,
                       typename Result::Polynomials{
                           /*sin=*/product(polynomials.sin, right.aperiodic_),
                           /*cos=*/product(polynomials.cos, right.aperiodic_)});
  }
  for (auto const& [ω, polynomials] : right.periodic_) {
    terms.emplace_back(ω,
                       typename Result::Polynomials{
                           /*sin=*/product(left.aperiodic_, polynomials.sin),
                           /*cos=*/product(left.aperiodic_, polynomials.cos)});
  }
  for (auto const& [ωl, polynomials_left] : left.periodic_) {
    for (auto const& [ωr, polynomials_right] : right.periodic_) {
      terms.emplace_back(
          ωl - ωr,
          typename Result::Polynomials{
              /*sin=*/(product(-polynomials_left.cos, polynomials_right.sin) +
                       product(polynomials_left.sin, polynomials_right.cos)) /
                  2,
              /*cos=*/(product(polynomials_left.sin, polynomials_right.sin) +
                       product(polynomials_left.cos, polynomials_right.cos)) /
                  2});
      terms.emplace_back(
          ωl + ωr,
          typename Result::Polynomials{
              /*sin=*/(product(polynomials_left.cos, polynomials_right.sin) +
                       product(polynomials_left.sin, polynomials_right.cos)) /
                  2,
              /*cos=*/(product(-polynomials_left.sin, polynomials_right.sin) +
                       product(polynomials_left.cos, polynomials_right.cos)) /
                  2});
    }
  }

  // Sort the vector by ascending frequency.
  std::sort(terms.begin(), terms.end(),
            [](Term const& left, Term const& right) {
              return left.first < right.first;
            });

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
      periodic.insert(periodic.cend(), *it);
    }
    previous_ω = ω;
  }

  return Result(std::move(aperiodic), std::move(periodic));
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
Instant const& PoissonSeries<Value, degree_, Evaluator>::origin() const {
  return origin_;
}

template<typename Value, int degree_,
         template<typename, typename, int> class Evaluator>
Value PoissonSeries<Value, degree_, Evaluator>::operator()(
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
PoissonSeries<Value, degree_, Evaluator>
PoissonSeries<Value, degree_, Evaluator>::AtOrigin(
    Instant const& origin) const {
  Time const shift = origin - origin_;
  auto aperiodic = aperiodic_.AtOrigin(origin);

  PolynomialsByAngularFrequency periodic;
  for (auto const& [ω, polynomials] : periodic_) {
    Polynomial const sin = polynomials.sin.AtOrigin(origin);
    Polynomial const cos = polynomials.cos.AtOrigin(origin);
    periodic.emplace_hint(
        periodic.cend(),
        ω,
        Polynomials{/*sin=*/sin * Cos(ω * shift) - cos * Sin(ω * shift),
                    /*cos=*/sin * Sin(ω * shift) + cos * Cos(ω * shift)});
  }
  return {std::move(aperiodic), std::move(periodic)};
}

template<typename Value, int degree_,
         template<typename, typename, int> class Evaluator>
PoissonSeries<quantities::Primitive<Value, Time>, degree_ + 1, Evaluator>
PoissonSeries<Value, degree_, Evaluator>::Primitive() const {
  using Result =
      PoissonSeries<quantities::Primitive<Value, Time>, degree_ + 1, Evaluator>;
  typename Result::PolynomialsByAngularFrequency periodic;
  typename Result::Polynomial aperiodic = aperiodic_.Primitive();
  for (auto const& [ω, polynomials] : periodic_) {
    periodic.emplace_hint(
        periodic.cend(),
        ω,
        AngularFrequencyPrimitive<Value, degree_, Evaluator>(ω, polynomials));
  }
  return Result{std::move(aperiodic), std::move(periodic)};
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
      aperiodic_primitive.Evaluate(t2) - aperiodic_primitive.Evaluate(t1);
  for (auto const& [ω, polynomials] : periodic_) {
    // Integration by parts.
    FirstPart const first_part(
        typename FirstPart::Polynomial({}, origin_),
        {{ω,
          {/*sin=*/typename FirstPart::Polynomial(polynomials.cos ),
           /*cos=*/typename FirstPart::Polynomial(-polynomials.sin)}}});
    result += (first_part(t2) - first_part(t1)) / ω * Radian;

    if constexpr (degree_ != 0) {
      auto const sin_polynomial =
          -polynomials.cos.template Derivative<1>();
      auto const cos_polynomial =
          polynomials.sin.template Derivative<1>();
      SecondPart const second_part(typename SecondPart::Polynomial({}, origin_),
                                   {{ω,
                                     {/*sin=*/sin_polynomial,
                                      /*cos=*/cos_polynomial}}});
      result += second_part.Integrate(t1, t2) / ω * Radian;
    }
  }
  return result;
}

template<typename Value, int degree_,
         template<typename, typename, int> class Evaluator>
template<typename V, int d, template<typename, typename, int> class E>
PoissonSeries<Value, degree_, Evaluator>&
PoissonSeries<Value, degree_, Evaluator>::operator+=(
    PoissonSeries<V, d, E> const& right) {
  *this = *this + right;
  return *this;
}

template<typename Value, int degree_,
         template<typename, typename, int> class Evaluator>
template<typename V, int d, template<typename, typename, int> class E>
PoissonSeries<Value, degree_, Evaluator>&
PoissonSeries<Value, degree_, Evaluator>::operator-=(
    PoissonSeries<V, d, E> const& right) {
  *this = *this - right;
  return *this;
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
  auto aperiodic = -right.aperiodic_;
  for (auto const& [ω, polynomials] : right.periodic_) {
    periodic.emplace_hint(
        periodic.cend(),
        ω,
        typename Result::Polynomials{/*sin=*/-polynomials.sin,
                                     /*cos=*/-polynomials.cos});
  }
  return {std::move(aperiodic), std::move(periodic)};
}

template<typename Value, int ldegree_, int rdegree_,
         template<typename, typename, int> class Evaluator>
PoissonSeries<Value, std::max(ldegree_, rdegree_), Evaluator>
operator+(PoissonSeries<Value, ldegree_, Evaluator> const& left,
          PoissonSeries<Value, rdegree_, Evaluator> const& right) {
  using Result = PoissonSeries<Value, std::max(ldegree_, rdegree_), Evaluator>;
  typename Result::PolynomialsByAngularFrequency periodic;
  auto aperiodic = left.aperiodic_ + right.aperiodic_;
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
      periodic.emplace_hint(
          periodic.cend(),
          ωl,
          typename Result::Polynomials{
              /*sin=*/typename Result::Polynomial(polynomials_left.sin),
              /*cos=*/typename Result::Polynomial(polynomials_left.cos)});
      ++it_left;
    } else if (ωr < ωl) {
      auto const& polynomials_right = it_right->second;
      periodic.emplace_hint(
          periodic.cend(),
          ωr,
          typename Result::Polynomials{
              /*sin=*/typename Result::Polynomial(polynomials_right.sin),
              /*cos=*/typename Result::Polynomial(polynomials_right.cos)});
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
  return {std::move(aperiodic), std::move(periodic)};
}

template<typename Value, int ldegree_, int rdegree_,
         template<typename, typename, int> class Evaluator>
PoissonSeries<Value, std::max(ldegree_, rdegree_), Evaluator>
operator-(PoissonSeries<Value, ldegree_, Evaluator> const& left,
          PoissonSeries<Value, rdegree_, Evaluator> const& right) {
  using Result = PoissonSeries<Value, std::max(ldegree_, rdegree_), Evaluator>;
  typename Result::PolynomialsByAngularFrequency periodic;
  auto aperiodic = left.aperiodic_ - right.aperiodic_;
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
      periodic.emplace_hint(
          periodic.cend(),
          ωl,
          typename Result::Polynomials{
              /*sin=*/typename Result::Polynomial(polynomials_left.sin),
              /*cos=*/typename Result::Polynomial(polynomials_left.cos)});
      ++it_left;
    } else if (ωr < ωl) {
      auto const& polynomials_right = it_right->second;
      periodic.emplace_hint(
          periodic.cend(),
          ωr,
          typename Result::Polynomials{
              /*sin=*/typename Result::Polynomial(-polynomials_right.sin),
              /*cos=*/typename Result::Polynomial(-polynomials_right.cos)});
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
  return {std::move(aperiodic), std::move(periodic)};
}

template<typename Scalar, typename Value, int degree_,
         template<typename, typename, int> class Evaluator>
PoissonSeries<Product<Scalar, Value>, degree_, Evaluator>
operator*(Scalar const& left,
          PoissonSeries<Value, degree_, Evaluator> const& right) {
  using Result = PoissonSeries<Product<Scalar, Value>, degree_, Evaluator>;
  typename Result::PolynomialsByAngularFrequency periodic;
  auto aperiodic = left * right.aperiodic_;
  for (auto const& [ω, polynomials] : right.periodic_) {
    periodic.emplace_hint(
        periodic.cend(),
        ω,
        typename Result::Polynomials{/*sin=*/left * polynomials.sin,
                                     /*cos=*/left * polynomials.cos});
  }
  return {std::move(aperiodic), std::move(periodic)};
}

template<typename Scalar, typename Value, int degree_,
         template<typename, typename, int> class Evaluator>
PoissonSeries<Product<Value, Scalar>, degree_, Evaluator>
operator*(PoissonSeries<Value, degree_, Evaluator> const& left,
          Scalar const& right) {
  using Result = PoissonSeries<Product<Scalar, Value>, degree_, Evaluator>;
  typename Result::PolynomialsByAngularFrequency periodic;
  auto aperiodic = left.aperiodic_ * right;
  for (auto const& [ω, polynomials] : left.periodic_) {
    periodic.emplace_hint(
        periodic.cend(),
        ω,
        typename Result::Polynomials{/*sin=*/polynomials.sin * right,
                                     /*cos=*/polynomials.cos * right});
  }
  return {std::move(aperiodic), std::move(periodic)};
}

template<typename Scalar, typename Value, int degree_,
         template<typename, typename, int> class Evaluator>
PoissonSeries<Quotient<Value, Scalar>, degree_, Evaluator>
operator/(PoissonSeries<Value, degree_, Evaluator> const& left,
          Scalar const& right) {
  using Result = PoissonSeries<Product<Scalar, Value>, degree_, Evaluator>;
  typename Result::PolynomialsByAngularFrequency periodic;
  auto aperiodic = left.aperiodic_ / right;
  for (auto const& [ω, polynomials] : left.periodic_) {
    periodic.emplace_hint(
        periodic.cend(),
        ω,
        typename Result::Polynomials{/*sin=*/polynomials.sin / right,
                                     /*cos=*/polynomials.cos / right});
  }
  return {std::move(aperiodic), std::move(periodic)};
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
typename Hilbert<LValue, RValue>::InnerProductType
Dot(PoissonSeries<LValue, ldegree_, Evaluator> const& left,
    PoissonSeries<RValue, rdegree_, Evaluator> const& right,
    PoissonSeries<double, wdegree_, Evaluator> const& weight,
    Instant const& t_min,
    Instant const& t_max) {
  auto const integrand = PointwiseInnerProduct(left, right) * weight;
  auto const primitive = integrand.Primitive();
  return (primitive(t_max) - primitive(t_min)) / (t_max - t_min);
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
Value PiecewisePoissonSeries<Value, degree_, Evaluator>::operator()(
    Instant const& t) const {
  if (t == bounds_.back()) {
    return series_.back()(t);
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
  return series_[it - bounds_.cbegin() - 1](t);
}

template<typename Value, int degree_,
         template<typename, typename, int> class Evaluator>
template<typename V, int d, template<typename, typename, int> class E>
PiecewisePoissonSeries<Value, degree_, Evaluator>&
PiecewisePoissonSeries<Value, degree_, Evaluator>::operator+=(
    PoissonSeries<V, d, E> const& right) {
  *this = *this + right;
  return *this;
}

template<typename Value, int degree_,
         template<typename, typename, int> class Evaluator>
template<typename V, int d, template<typename, typename, int> class E>
PiecewisePoissonSeries<Value, degree_, Evaluator>&
PiecewisePoissonSeries<Value, degree_, Evaluator>::operator-=(
    PoissonSeries<V, d, E> const& right) {
  *this = *this - right;
  return *this;
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

// In practice changing the origin of the piecewise series chunks is horribly
// ill-conditioned, so the code below changes the origin of the (single) Poisson
// series.
// TODO(phl): All these origin changes might be expensive, see if we can factor
// them.

template<typename Value, int ldegree_, int rdegree_,
         template<typename, typename, int> class Evaluator>
PiecewisePoissonSeries<Value, std::max(ldegree_, rdegree_), Evaluator>
operator+(PoissonSeries<Value, ldegree_, Evaluator> const& left,
          PiecewisePoissonSeries<Value, rdegree_, Evaluator> const& right) {
  using Result =
      PiecewisePoissonSeries<Value, std::max(ldegree_, rdegree_), Evaluator>;
  std::vector<typename Result::Series> series;
  for (int i = 0; i < right.series_.size(); ++i) {
    Instant const origin = right.series_[i].origin();
    series.push_back(left.AtOrigin(origin) + right.series_[i]);
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
    Instant const origin = left.series_[i].origin();
    series.push_back(left.series_[i] + right.AtOrigin(origin));
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
    Instant const origin = right.series_[i].origin();
    series.push_back(left.AtOrigin(origin) - right.series_[i]);
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
    Instant const origin = left.series_[i].origin();
    series.push_back(left.series_[i] - right.AtOrigin(origin));
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
    Instant const origin = right.series_[i].origin();
    series.push_back(left.AtOrigin(origin) * right.series_[i]);
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
    Instant const origin = left.series_[i].origin();
    series.push_back(left.series_[i] * right.AtOrigin(origin));
  }
  return Result(left.bounds_, series);
}

template<typename LValue, typename RValue,
         int ldegree_, int rdegree_, int wdegree_,
         template<typename, typename, int> class Evaluator>
typename Hilbert<LValue, RValue>::InnerProductType
Dot(PoissonSeries<LValue, ldegree_, Evaluator> const& left,
    PiecewisePoissonSeries<RValue, rdegree_, Evaluator> const& right,
    PoissonSeries<double, wdegree_, Evaluator> const& weight) {
  return Dot(left, right, weight, right.t_min(), right.t_max());
}

template<typename LValue, typename RValue,
         int ldegree_, int rdegree_, int wdegree_,
         template<typename, typename, int> class Evaluator>
typename Hilbert<LValue, RValue>::InnerProductType
Dot(PoissonSeries<LValue, ldegree_, Evaluator> const& left,
    PiecewisePoissonSeries<RValue, rdegree_, Evaluator> const& right,
    PoissonSeries<double, wdegree_, Evaluator> const& weight,
    Instant const& t_min,
    Instant const& t_max) {
  using Result =
      Primitive<typename Hilbert<LValue, RValue>::InnerProductType, Time>;
  Result result{};
  for (int i = 0; i < right.series_.size(); ++i) {
    Instant const origin = right.series_[i].origin();
    auto const integrand =
        PointwiseInnerProduct(left.AtOrigin(origin), right.series_[i]) *
        weight.AtOrigin(origin);
    auto const integral =
        quadrature::GaussLegendre<ldegree_ + rdegree_ + wdegree_>(
            integrand, right.bounds_[i], right.bounds_[i + 1]);
    result += integral;
  }
  return result / (t_max - t_min);
}

template<typename LValue, typename RValue,
         int ldegree_, int rdegree_, int wdegree_,
         template<typename, typename, int> class Evaluator>
typename Hilbert<LValue, RValue>::InnerProductType
Dot(PiecewisePoissonSeries<LValue, ldegree_, Evaluator> const& left,
    PoissonSeries<RValue, rdegree_, Evaluator> const& right,
    PoissonSeries<double, wdegree_, Evaluator> const& weight) {
  return Dot(left, right, weight, left.t_min(), left.t_max());
}

template<typename LValue, typename RValue,
         int ldegree_, int rdegree_, int wdegree_,
         template<typename, typename, int> class Evaluator>
typename Hilbert<LValue, RValue>::InnerProductType
Dot(PiecewisePoissonSeries<LValue, ldegree_, Evaluator> const& left,
    PoissonSeries<RValue, rdegree_, Evaluator> const& right,
    PoissonSeries<double, wdegree_, Evaluator> const& weight,
    Instant const& t_min,
    Instant const& t_max) {
  using Result =
      Primitive<typename Hilbert<LValue, RValue>::InnerProductType, Time>;
  Result result{};
#if 0
  if (do_the_logging) {
      frequency_analysis::logger.Append(
          "right",
          right,
          mathematica::ExpressIn(Metre, Radian, Second));
      frequency_analysis::logger.Append(
          "weight",
          weight,
          mathematica::ExpressIn(Metre, Radian, Second));
  }
#endif
  for (int i = 0; i < left.series_.size(); ++i) {
    Instant const origin = left.series_[i].origin();
    auto const integrand =
        PointwiseInnerProduct(left.series_[i], right.AtOrigin(origin)) *
        weight.AtOrigin(origin);
    auto const integral =
        quadrature::GaussLegendre<ldegree_ + rdegree_ + wdegree_>(
            integrand, left.bounds_[i], left.bounds_[i + 1]);
#if 0
    if (do_the_logging) {
      frequency_analysis::logger.Append(
          "rightAtOrigin",
          right.AtOrigin(origin),
          mathematica::ExpressIn(Metre, Radian, Second));
      frequency_analysis::logger.Append(
          "weigthAtOrigin",
          weight.AtOrigin(origin),
          mathematica::ExpressIn(Metre, Radian, Second));
      frequency_analysis::logger.Append(
          "series",
          left.series_[i],
          mathematica::ExpressIn(Metre, Radian, Second));
      frequency_analysis::logger.Append(
          "integrand",
          integrand,
          mathematica::ExpressIn(Metre, Radian, Second));
      frequency_analysis::logger.Append(
          "bounds",
          std::vector{left.bounds_[i + 1], left.bounds_[i]},
          mathematica::ExpressIn(Metre, Radian, Second));
      frequency_analysis::logger.Append(
          "primitive",
          primitive,
          mathematica::ExpressIn(Metre, Radian, Second));
      frequency_analysis::logger.Append(
          "integral",
          integral.value,
          mathematica::ExpressIn(Metre, Radian, Second));
    }
#endif
    result += integral;
  }
  return result / (t_max - t_min);
}

}  // namespace internal_poisson_series
}  // namespace numerics
}  // namespace principia
