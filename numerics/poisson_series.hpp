
#pragma once

#include <algorithm>
#include <map>

#include "numerics/polynomial.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"

namespace principia {
namespace numerics {
namespace internal_poisson_series {

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
  using Polynomial =
      numerics::PolynomialInMonomialBasis<Value, Time, degree_, Evaluator>;

  // TODO(phl): Use designated initializers for this struct once this project
  // can be compiled using c++latest.
  struct Polynomials {
    Polynomial sin;
    Polynomial cos;
  };

  using PolynomialsByAngularFrequency = std::map<AngularFrequency, Polynomials>;

  PoissonSeries(Polynomial const& aperiodic,
                PolynomialsByAngularFrequency const& periodic);

  Value Evaluate(Time const& t) const;

  // The constant term of the result is zero.
  PoissonSeries<Primitive<Value, Time>, degree_ + 1, Evaluator> Primitive()
      const;

 private:
  static Polynomials AngularFrequencyPrimitive(AngularFrequency const& ω,
                                               Polynomials const& polynomials);

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

}  // namespace internal_poisson_series

using internal_poisson_series::PoissonSeries;

}  // namespace numerics
}  // namespace principia

#include "numerics/poisson_series_body.hpp"
