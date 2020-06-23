
#pragma once

#include "numerics/poisson_series.hpp"

#include <functional>

#include "quantities/elementary_functions.hpp"

namespace principia {
namespace numerics {
namespace internal_poisson_series {

using quantities::Cos;
using quantities::Infinity;
using quantities::Sin;

template<typename Value, int degree_,
         template<typename, typename, int> class Evaluator>
PoissonSeries<Value, degree_, Evaluator>::PoissonSeries(
    Polynomial const& aperiodic,
    PolynomialsByAngularFrequency const& periodic)
    : aperiodic_(aperiodic), periodic_(periodic) {
  for (auto const& [ω, _] : periodic_) {
    CHECK_LT(AngularFrequency{}, ω);
  }
}

template<typename Value, int degree_,
         template<typename, typename, int> class Evaluator>
Value PoissonSeries<Value, degree_, Evaluator>::Evaluate(Time const& t) const {
  Value result = aperiodic_.Evaluate(t);
  for (auto const& [ω, polynomials] : periodic_) {
    result += polynomials.sin.Evaluate(t) * Sin(ω * t) +
              polynomials.cos.Evaluate(t) * Cos(ω * t);
  }
  return result;
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
  PoissonSeries<Value, rdegree_, Evaluator> result;
  result.aperiodic_ = -right.aperiodic_;
  for (auto const& [ω, polynomials] : right.periodic_) {
    result.emplace(ω, {.sin = -polynomials.sin, .cos = -polynomials.cos});
  }
  return result;
}

template<typename Value, int ldegree_, int rdegree_,
         template<typename, typename, int> class Evaluator>
PoissonSeries<Value, std::max(ldegree_, rdegree_), Evaluator>
operator+(PoissonSeries<Value, ldegree_, Evaluator> const& left,
          PoissonSeries<Value, rdegree_, Evaluator> const& right) {
  PoissonSeries<Value, std::max(ldegree_, rdegree_), Evaluator> result;
  result.aperiodic_ = left.aperiodic_ + right.aperiodic_;
  auto it_left = left.periodic_.cbegin();
  auto it_right = right.periodic_.cbegin();
  while (it_left != left.periodic_.cend() && it_right != right.periodic_.cend()) {
    auto const ωl = it_left.first;
    auto const ωr = it_right.first;
    if (it_right == right.periodic_.cend() || ωl < ωr) {
      result.periodic_.insert(*it_left);
      ++it_left;
    } else if (it_left == left.periodic_.cend() || ωr < ωl) {
      result.periodic_.insert(*it_right);
      ++it_right;
    } else {
      DCHECK_EQ(ωl, ωr);
      auto const& polynomials_left = it_left->second;
      auto const& polynomials_right = it_right->second;
      result.periodic_.emplace(
          ωl,
          {.sin = polynomials_left.sin + polynomials_right.sin,
           .cos = polynomials_left.con + polynomials_right.cos});
    }
  }
  return result;
}

template<typename Value, int ldegree_, int rdegree_,
         template<typename, typename, int> class Evaluator>
PoissonSeries<Value, std::max(ldegree_, rdegree_), Evaluator>
operator-(PoissonSeries<Value, ldegree_, Evaluator> const& left,
          PoissonSeries<Value, rdegree_, Evaluator> const& right) {
  PoissonSeries<Value, std::max(ldegree_, rdegree_), Evaluator> result;
  result.aperiodic_ = left.aperiodic_ + right.aperiodic_;
  auto it_left = left.periodic_.cbegin();
  auto it_right = right.periodic_.cbegin();
  while (it_left != left.periodic_.cend() && it_right != right.periodic_.cend()) {
    auto const ωl = it_left.first;
    auto const ωr = it_right.first;
    if (it_right == right.periodic_.cend() || ωl < ωr) {
      result.periodic_.insert(*it_left);
      ++it_left;
    } else if (it_left == left.periodic_.cend() || ωr < ωl) {
      auto const& polynomials_right = it_right->second;
      result.periodic_.emplace(ωr,
                               {.sin = -polynomials_right.sin,
                                .cos = -polynomials_right.cos});
      ++it_right;
    } else {
      DCHECK_EQ(ωl, ωr);
      auto const& polynomials_left = it_left->second;
      auto const& polynomials_right = it_right->second;
      result.periodic_.emplace(
          ωl,
          {.sin = polynomials_left.sin - polynomials_right.sin,
           .cos = polynomials_left.con - polynomials_right.cos});
    }
  }
  return result;
}

template<typename Scalar, typename Value, int degree_,
         template<typename, typename, int> class Evaluator>
PoissonSeries<Product<Scalar, Value>, degree_, Evaluator>
operator*(Scalar const& left,
          PoissonSeries<Value, degree_, Evaluator> const& right) {
  PoissonSeries<Value, rdegree_, Evaluator> result;
  result.aperiodic_ = left * right.aperiodic_;
  for (auto const& [ω, polynomials] : right.periodic_) {
    result.emplace(ω,
                   {.sin = left * polynomials.sin,
                    .cos = left * polynomials.cos});
  }
  return result;
}

template<typename Scalar, typename Value, int degree_,
         template<typename, typename, int> class Evaluator>
PoissonSeries<Product<Value, Scalar>, degree_, Evaluator>
operator*(PoissonSeries<Value, degree_, Evaluator> const& left,
          Scalar const& right) {
  PoissonSeries<Value, rdegree_, Evaluator> result;
  result.aperiodic_ = left.aperiodic_ * right;
  for (auto const& [ω, polynomials] : left.periodic_) {
    result.emplace(ω,
                   {.sin = polynomials.sin * right,
                    .cos = polynomials.cos * right});
  }
  return result;
}

template<typename Scalar, typename Value, int degree_,
         template<typename, typename, int> class Evaluator>
PoissonSeries<Quotient<Value, Scalar>, degree_, Evaluator>
operator/(PoissonSeries<Value, degree_, Evaluator> const& left,
          Scalar const& right) {
  PoissonSeries<Value, rdegree_, Evaluator> result;
  result.aperiodic_ = left.aperiodic_ / right;
  for (auto const& [ω, polynomials] : left.periodic_) {
    result.emplace(ω,
                   {.sin = polynomials.sin / right,
                    .cos = polynomials.cos / right});
  }
  return result;
}

}  // namespace internal_poisson_series
}  // namespace numerics
}  // namespace principia
