
#pragma once

#include "numerics/poisson_series.hpp"

#include <algorithm>
#include <functional>

#include "quantities/elementary_functions.hpp"

namespace principia {
namespace numerics {
namespace internal_poisson_series {

using quantities::Cos;
using quantities::Infinity;
using quantities::Sin;

// TODO(phl): This file should probably take advantage of emplace_hint.

template<typename Value, int degree_,
         template<typename, typename, int> class Evaluator>
PoissonSeries<Value, degree_, Evaluator>::PoissonSeries(
    Polynomial const& aperiodic,
    PolynomialsByAngularFrequency const& periodic)
    : aperiodic_(aperiodic) {
  // The |periodic| map may have elements with positive or negative angular
  // frequencies.  Normalize our member variable to only have positive angular
  // frequencies.
  // TODO(phl): It would be good to have +=, -=, etc. for polynomials.
  for (auto it = periodic.crbegin(); it != periodic.crend(); ++it) {
    auto const ω = it->first;
    if (ω < AngularFrequency{}) {
      auto const positive_it = periodic_.find(-ω);
      if (positive_it == periodic_.cend()) {
        periodic_.emplace(-ω,
                          Polynomials{/*sin=*/-it->second.sin,
                                      /*cos=*/it->second.cos});
      } else {
        positive_it->second = {
            /*sin=*/positive_it->second.sin - it->second.sin,
            /*cos=*/positive_it->second.cos + it->second.cos};
      }
    } else if (ω > AngularFrequency{}) {
      periodic_.insert(*it);
    } else {
      aperiodic_ = aperiodic_ + it->second.cos;
    }
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
  using Result = PoissonSeries<Value, rdegree_, Evaluator>;
  typename Result::PolynomialsByAngularFrequency periodic;
  auto const aperiodic = -right.aperiodic_;
  for (auto const& [ω, polynomials] : right.periodic_) {
    periodic.emplace(ω,
                     Result::Polynomials{/*sin=*/-polynomials.sin,
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
      periodic.insert(*it_left);
      ++it_left;
    } else if (ωr < ωl) {
      periodic.insert(*it_right);
      ++it_right;
    } else {
      DCHECK_EQ(ωl, ωr);
      auto const& polynomials_left = it_left->second;
      auto const& polynomials_right = it_right->second;
      periodic.emplace(
          ωl,
          Result::Polynomials{
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
      periodic.insert(*it_left);
      ++it_left;
    } else if (ωr < ωl) {
      auto const& polynomials_right = it_right->second;
      periodic.emplace(ωr,
                       Result::Polynomials{/*sin=*/-polynomials_right.sin,
                                           /*cos=*/-polynomials_right.cos});
      ++it_right;
    } else {
      DCHECK_EQ(ωl, ωr);
      auto const& polynomials_left = it_left->second;
      auto const& polynomials_right = it_right->second;
      periodic.emplace(
          ωl,
          Result::Polynomials{
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
    periodic.emplace(ω,
                     Result::Polynomials{/*sin=*/left * polynomials.sin,
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
    periodic.emplace(ω,
                     Result::Polynomials{/*sin=*/polynomials.sin * right,
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
    periodic.emplace(ω,
                     Result::Polynomials{/*sin=*/polynomials.sin / right,
                                         /*cos=*/polynomials.cos / right});
  }
  return {aperiodic, std::move(periodic)};
}

}  // namespace internal_poisson_series
}  // namespace numerics
}  // namespace principia
