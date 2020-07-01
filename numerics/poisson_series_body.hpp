
#pragma once

#include "numerics/poisson_series.hpp"

#include <algorithm>
#include <functional>
#include <map>
#include <optional>

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
    : aperiodic_(aperiodic) {
  // The |periodic| map may have elements with positive or negative angular
  // frequencies.  Normalize our member variable to only have positive angular
  // frequencies.
  for (auto it = periodic.crbegin(); it != periodic.crend(); ++it) {
    auto const ω = it->first;
    if (ω < AngularFrequency{}) {
      auto const positive_it = periodic_.find(-ω);
      if (positive_it == periodic_.cend()) {
        periodic_.emplace(-ω,
                          Polynomials{/*sin=*/-it->second.sin,
                                      /*cos=*/it->second.cos});
      } else {
        positive_it->second.sin -= it->second.sin;
        positive_it->second.cos += it->second.cos;
      }
    } else if (ω > AngularFrequency{}) {
      periodic_.insert(periodic_.cbegin(), *it);
    } else {
      aperiodic_ += it->second.cos;
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

template<typename Value, int degree_,
         template<typename, typename, int> class Evaluator>
PoissonSeries<Primitive<Value, Time>, degree_ + 1, Evaluator>
PoissonSeries<Value, degree_, Evaluator>::Primitive() const {
  using Result =
      PoissonSeries<quantities::Primitive<Value, Time>, degree_ + 1, Evaluator>;
  typename Result::PolynomialsByAngularFrequency periodic;
  typename Result::Polynomial const aperiodic = aperiodic_.Primitive();
  for (auto const& [ω, polynomials] : periodic_) {
    periodic.emplace_hint(
        periodic.cend(),
        ω,
        AngularFrequencyPrimitive(ω, polynomials));
  }
  return Result{aperiodic, std::move(periodic)};
}

template<typename Value, int degree_,
         template<typename, typename, int> class Evaluator>
typename PoissonSeries<Value, degree_, Evaluator>::Polynomials
PoissonSeries<Value, degree_, Evaluator>::AngularFrequencyPrimitive(
    AngularFrequency const& ω,
    Polynomials const& polynomials) {
  // Integration by parts.
  auto const first_part = Polynomials{/*sin=*/polynomials.cos / ω,
                                      /*cos=*/-polynomials.sin / ω};
  if constexpr (degree_ == 0) {
    return first_part;
  } else {
    auto const second_part = AngularFrequencyPrimitive(
        ω,
        Polynomials{/*sin=*/-polynomials.cos.Derivative<1>() / ω,
                    /*cos=*/polynomials.sin.Derivative<1>() / ω});
    return Polynomials{/*sin=*/first_part.sin + second_part.sin,
                       /*cos=*/first_part.cos + second_part.cos};
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

}  // namespace internal_poisson_series
}  // namespace numerics
}  // namespace principia
