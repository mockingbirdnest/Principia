
#pragma once

#include <algorithm>
#include <map>
#include <vector>

#include "geometry/interval.hpp"
#include "geometry/named_quantities.hpp"
#include "numerics/polynomial.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"

namespace principia {
namespace numerics {
namespace internal_poisson_series {

using geometry::Instant;
using geometry::Interval;
using quantities::AngularFrequency;
using quantities::Primitive;
using quantities::Product;
using quantities::Quotient;
using quantities::Time;

// A Poisson series is the sum of terms of the form:
//   aₙtⁿ      aₙₖ tⁿ sin ωₖ t      aₙₖ tⁿ cos ωₖ t
// Terms of the first kind are called aperiodic, terms of the second and third
// kind are called periodic.  Poisson series form an algebra that is stable by
// derivation and integration.
template<typename Value, int degree_,
         template<typename, typename, int> class Evaluator>
class PoissonSeries {
 public:
  static const int degree = degree_;
  using Polynomial =
      numerics::PolynomialInMonomialBasis<Value, Instant, degree_, Evaluator>;

  // TODO(phl): Use designated initializers for this struct once this project
  // can be compiled using c++latest.
  struct Polynomials {
    Polynomial sin;
    Polynomial cos;
  };

  using PolynomialsByAngularFrequency = std::map<AngularFrequency, Polynomials>;

  PoissonSeries(Polynomial const& aperiodic,
                PolynomialsByAngularFrequency const& periodic);

  Instant const& origin() const;

  Value operator()(Instant const& t) const;

  // The constant term of the result is zero.
  PoissonSeries<quantities::Primitive<Value, Time>, degree_ + 1, Evaluator>
  Primitive() const;

 private:
  Instant origin_;  // Common to all polynomials.
  Polynomial aperiodic_;
  // All the keys in this map are positive.
  PolynomialsByAngularFrequency periodic_;

  template<typename V, int r, template<typename, typename, int> class E>
  PoissonSeries<V, r, E>
  friend operator-(PoissonSeries<V, r, E> const& right);
  template<typename V, int l, int r, template<typename, typename, int> class E>
  PoissonSeries<V, std::max(l, r), E>
  friend operator+(PoissonSeries<V, l, E> const& left,
                   PoissonSeries<V, r, E> const& right);
  template<typename V, int l, int r, template<typename, typename, int> class E>
  PoissonSeries<V, std::max(l, r), E>
  friend operator-(PoissonSeries<V, l, E> const& left,
                   PoissonSeries<V, r, E> const& right);
  template<typename Scalar,
           typename V, int d,
           template<typename, typename, int> class E>
  PoissonSeries<Product<Scalar, V>, d, E>
  friend operator*(Scalar const& left,
                   PoissonSeries<V, d, E> const& right);
  template<typename Scalar,
           typename V, int d,
           template<typename, typename, int> class E>
  PoissonSeries<Product<V, Scalar>, d, E>
  friend operator*(PoissonSeries<V, d, E> const& left,
                   Scalar const& right);
  template<typename Scalar,
           typename V, int d,
           template<typename, typename, int> class E>
  PoissonSeries<Quotient<V, Scalar>, d, E>
  friend operator/(PoissonSeries<V, d, E> const& left,
                   Scalar const& right);
  template<typename L, typename R,
           int l, int r,
           template<typename, typename, int> class E>
  PoissonSeries<Product<L, R>, l + r, E>
  friend operator*(PoissonSeries<L, l, E> const& left,
                   PoissonSeries<R, r, E> const& right);
};

// Vector space of Poisson series.

template<typename Value, int rdegree_,
         template<typename, typename, int> class Evaluator>
PoissonSeries<Value, rdegree_, Evaluator>
operator+(PoissonSeries<Value, rdegree_, Evaluator> const& right);

template<typename Value, int rdegree_,
         template<typename, typename, int> class Evaluator>
PoissonSeries<Value, rdegree_, Evaluator>
operator-(PoissonSeries<Value, rdegree_, Evaluator> const& right);

template<typename Value, int ldegree_, int rdegree_,
         template<typename, typename, int> class Evaluator>
PoissonSeries<Value, std::max(ldegree_, rdegree_), Evaluator>
operator+(PoissonSeries<Value, ldegree_, Evaluator> const& left,
          PoissonSeries<Value, rdegree_, Evaluator> const& right);

template<typename Value, int ldegree_, int rdegree_,
         template<typename, typename, int> class Evaluator>
PoissonSeries<Value, std::max(ldegree_, rdegree_), Evaluator>
operator-(PoissonSeries<Value, ldegree_, Evaluator> const& left,
          PoissonSeries<Value, rdegree_, Evaluator> const& right);

template<typename Scalar,
         typename Value, int degree_,
         template<typename, typename, int> class Evaluator>
PoissonSeries<Product<Scalar, Value>, degree_, Evaluator>
operator*(Scalar const& left,
          PoissonSeries<Value, degree_, Evaluator> const& right);

template<typename Scalar,
         typename Value, int degree_,
         template<typename, typename, int> class Evaluator>
PoissonSeries<Product<Value, Scalar>, degree_, Evaluator>
operator*(PoissonSeries<Value, degree_, Evaluator> const& left,
          Scalar const& right);

template<typename Scalar,
         typename Value, int degree_,
         template<typename, typename, int> class Evaluator>
PoissonSeries<Quotient<Value, Scalar>, degree_, Evaluator>
operator/(PoissonSeries<Value, degree_, Evaluator> const& left,
          Scalar const& right);

// Algebra of Poisson series.

template<typename LValue, typename RValue,
         int ldegree_, int rdegree_,
         template<typename, typename, int> class Evaluator>
PoissonSeries<Product<LValue, RValue>, ldegree_ + rdegree_, Evaluator>
operator*(PoissonSeries<LValue, ldegree_, Evaluator> const& left,
          PoissonSeries<RValue, rdegree_, Evaluator> const& right);

// Inner product space of Poisson series.

// Technically the weight function must be nonnegative for this to be an inner
// product.  Not sure how this works with the flat-top windows, which can be
// negative.
template<typename LValue, typename RValue,
         int ldegree_, int rdegree_, int wdegree_,
         template<typename, typename, int> class Evaluator>
Primitive<Product<LValue, RValue>, Time>
Dot(PoissonSeries<LValue, ldegree_, Evaluator> const& left,
    PoissonSeries<RValue, rdegree_, Evaluator> const& right,
    PoissonSeries<double, wdegree_, Evaluator> const& weight,
    Instant const& t_min,
    Instant const& t_max);

// A function defined by Poisson series piecewise.  Each of the Poisson series
// making up the function applies over the semi-open interval
// [internal.min, interval.max[.  It's not required that the function be
// continuous.
template<typename Value, int degree_,
         template<typename, typename, int> class Evaluator>
class PiecewisePoissonSeries {
 public:
  using Series = PoissonSeries<Value, degree_, Evaluator>;

  PiecewisePoissonSeries(Interval<Instant> const& interval,
                         Series const& series);

  // The intervals for successive calls to Append must be consecutive.  For the
  // first call, the interval must be consecutive with the one passed at
  // construction.
  void Append(Interval<Instant> const& interval,
              Series const& series);

  Instant t_min() const;
  Instant t_max() const;

  // t must be in the interval [t_min, t_max[.
  Value operator()(Instant const& t) const;

 private:
  PiecewisePoissonSeries(std::vector<Instant> const& bounds,
                         std::vector<Series> const& series);

  std::vector<Instant> bounds_;
  std::vector<Series> series_;

  template<typename V, int r, template<typename, typename, int> class E>
  PiecewisePoissonSeries<V, r, E>
  friend operator-(PiecewisePoissonSeries<V, r, E> const& right);
  template<typename S, typename V, int d,
           template<typename, typename, int> class E>
  PiecewisePoissonSeries<Product<S, V>, d, E>
  friend operator*(S const& left,
                   PiecewisePoissonSeries<V, d, E> const& right);
  template<typename S, typename V, int d,
           template<typename, typename, int> class E>
  PiecewisePoissonSeries<Product<V, S>, d, E>
  friend operator*(PiecewisePoissonSeries<V, d, E> const& left,
                   S const& right);
  template<typename S, typename V, int d,
           template<typename, typename, int> class E>
  PiecewisePoissonSeries<Quotient<V, S>, d, E>
  friend operator/(PiecewisePoissonSeries<V, d, E> const& left,
                   S const& right);
  template<typename V, int l, int r, template<typename, typename, int> class E>
  PiecewisePoissonSeries<V, std::max(l, r), E>
  friend operator+(PoissonSeries<V, l, E> const& left,
                   PiecewisePoissonSeries<V, r, E> const& right);
  template<typename V, int l, int r, template<typename, typename, int> class E>
  PiecewisePoissonSeries<V, std::max(l, r), E>
  friend operator+(PiecewisePoissonSeries<V, l, E> const& left,
                   PoissonSeries<V, r, E> const& right);
  template<typename V, int l, int r, template<typename, typename, int> class E>
  PiecewisePoissonSeries<V, std::max(l, r), E>
  friend operator-(PoissonSeries<V, l, E> const& left,
                   PiecewisePoissonSeries<V, r, E> const& right);
  template<typename V, int l, int r, template<typename, typename, int> class E>
  PiecewisePoissonSeries<V, std::max(l, r), E>
  friend operator-(PiecewisePoissonSeries<V, l, E> const& left,
                   PoissonSeries<V, r, E> const& right);
  template<typename L, typename R, int l, int r,
           template<typename, typename, int> class E>
  PiecewisePoissonSeries<Product<L, R>, l + r, E>
  friend operator*(PoissonSeries<L, l, E> const& left,
                   PiecewisePoissonSeries<R, r, E> const& right);
  template<typename L, typename R, int l, int r,
           template<typename, typename, int> class E>
  PiecewisePoissonSeries<Product<L, R>, l + r, E>
  friend operator*(PiecewisePoissonSeries<L, l, E> const& left,
                   PoissonSeries<R, r, E> const& right);
  template<typename L, typename R, int l, int r, int w,
           template<typename, typename, int> class E>
  quantities::Primitive<Product<L, R>, Time>
  friend Dot(PoissonSeries<L, l, E> const& left,
             PiecewisePoissonSeries<R, r, E> const& right,
             PoissonSeries<double, w, E> const& weight);
  template<typename L, typename R, int l, int r, int w,
           template<typename, typename, int> class E>
  quantities::Primitive<Product<L, R>, Time>
  friend Dot(PiecewisePoissonSeries<L, l, E> const& left,
             PoissonSeries<R, r, E> const& right,
             PoissonSeries<double, w, E> const& weight);
};

// Some of the vector space operations for piecewise Poisson series.  Note that
// while piecewise Poisson series have an algebra structure, we do not need it.

template<typename Value, int rdegree_,
         template<typename, typename, int> class Evaluator>
PiecewisePoissonSeries<Value, rdegree_, Evaluator>
operator+(PiecewisePoissonSeries<Value, rdegree_, Evaluator> const& right);

template<typename Value, int rdegree_,
         template<typename, typename, int> class Evaluator>
PiecewisePoissonSeries<Value, rdegree_, Evaluator>
operator-(PiecewisePoissonSeries<Value, rdegree_, Evaluator> const& right);

template<typename Scalar,
         typename Value, int degree_,
         template<typename, typename, int> class Evaluator>
PiecewisePoissonSeries<Product<Scalar, Value>, degree_, Evaluator>
operator*(Scalar const& left,
          PiecewisePoissonSeries<Value, degree_, Evaluator> const& right);

template<typename Scalar,
         typename Value, int degree_,
         template<typename, typename, int> class Evaluator>
PiecewisePoissonSeries<Product<Value, Scalar>, degree_, Evaluator>
operator*(PiecewisePoissonSeries<Value, degree_, Evaluator> const& left,
          Scalar const& right);

template<typename Scalar,
         typename Value, int degree_,
         template<typename, typename, int> class Evaluator>
PiecewisePoissonSeries<Quotient<Value, Scalar>, degree_, Evaluator>
operator/(PiecewisePoissonSeries<Value, degree_, Evaluator> const& left,
          Scalar const& right);

// Action of Poisson series on piecewise Poisson series.

template<typename Value, int ldegree_, int rdegree_,
         template<typename, typename, int> class Evaluator>
PiecewisePoissonSeries<Value, std::max(ldegree_, rdegree_), Evaluator>
operator+(PoissonSeries<Value, ldegree_, Evaluator> const& left,
          PiecewisePoissonSeries<Value, rdegree_, Evaluator> const& right);

template<typename Value, int ldegree_, int rdegree_,
         template<typename, typename, int> class Evaluator>
PiecewisePoissonSeries<Value, std::max(ldegree_, rdegree_), Evaluator>
operator+(PiecewisePoissonSeries<Value, ldegree_, Evaluator> const& left,
          PoissonSeries<Value, rdegree_, Evaluator> const& right);

template<typename Value, int ldegree_, int rdegree_,
         template<typename, typename, int> class Evaluator>
PiecewisePoissonSeries<Value, std::max(ldegree_, rdegree_), Evaluator>
operator-(PoissonSeries<Value, ldegree_, Evaluator> const& left,
          PiecewisePoissonSeries<Value, rdegree_, Evaluator> const& right);

template<typename Value, int ldegree_, int rdegree_,
         template<typename, typename, int> class Evaluator>
PiecewisePoissonSeries<Value, std::max(ldegree_, rdegree_), Evaluator>
operator-(PiecewisePoissonSeries<Value, ldegree_, Evaluator> const& left,
          PoissonSeries<Value, rdegree_, Evaluator> const& right);

template<typename LValue, typename RValue,
         int ldegree_, int rdegree_,
         template<typename, typename, int> class Evaluator>
PiecewisePoissonSeries<Product<LValue, RValue>, ldegree_ + rdegree_, Evaluator>
operator*(PoissonSeries<LValue, ldegree_, Evaluator> const& left,
          PiecewisePoissonSeries<RValue, rdegree_, Evaluator> const& right);

template<typename LValue, typename RValue,
         int ldegree_, int rdegree_,
         template<typename, typename, int> class Evaluator>
PiecewisePoissonSeries<Product<LValue, RValue>, ldegree_ + rdegree_, Evaluator>
operator*(PiecewisePoissonSeries<LValue, ldegree_, Evaluator> const& left,
          PoissonSeries<RValue, rdegree_, Evaluator> const& right);

template<typename LValue, typename RValue,
         int ldegree_, int rdegree_, int wdegree_,
         template<typename, typename, int> class Evaluator>
Primitive<Product<LValue, RValue>, Time>
Dot(PoissonSeries<LValue, ldegree_, Evaluator> const& left,
    PiecewisePoissonSeries<RValue, rdegree_, Evaluator> const& right,
    PoissonSeries<double, wdegree_, Evaluator> const& weight);

template<typename LValue, typename RValue,
         int ldegree_, int rdegree_, int wdegree_,
         template<typename, typename, int> class Evaluator>
Primitive<Product<LValue, RValue>, Time>
Dot(PiecewisePoissonSeries<LValue, ldegree_, Evaluator> const& left,
    PoissonSeries<RValue, rdegree_, Evaluator> const& right,
    PoissonSeries<double, wdegree_, Evaluator> const& weight);

}  // namespace internal_poisson_series

using internal_poisson_series::PiecewisePoissonSeries;
using internal_poisson_series::PoissonSeries;

}  // namespace numerics
}  // namespace principia

#include "numerics/poisson_series_body.hpp"
