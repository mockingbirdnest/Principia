
#pragma once

#include <algorithm>
#include <functional>
#include <map>
#include <ostream>
#include <string>
#include <utility>
#include <vector>

#include "base/macros.hpp"
#include "base/not_null.hpp"
#include "geometry/complexification.hpp"
#include "geometry/hilbert.hpp"
#include "geometry/interval.hpp"
#include "geometry/named_quantities.hpp"
#include "numerics/polynomial.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"
#include "serialization/numerics.pb.h"

namespace principia {
namespace numerics {
FORWARD_DECLARE_FROM(poisson_series,
                     TEMPLATE(typename Value,
                              int aperiodic_degree, int periodic_degree,
                              template<typename, typename, int> class Evaluator)
                              class,
                     PoissonSeries);
}  // namespace numerics

namespace mathematica {
FORWARD_DECLARE_FUNCTION_FROM(
    mathematica,
    TEMPLATE(typename Value,
             int aperiodic_degree, int periodic_degree,
             template<typename, typename, int> class Evaluator,
             typename OptionalExpressIn) std::string,
    ToMathematicaExpression,
    (numerics::PoissonSeries<Value,
                             aperiodic_degree, periodic_degree,
                             Evaluator> const& series,
     OptionalExpressIn express_in));
}  // namespace mathematica

namespace numerics {
namespace internal_poisson_series {

using base::not_null;
using geometry::Complexification;
using geometry::Hilbert;
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
template<typename Value,
         int aperiodic_degree_, int periodic_degree_,
         template<typename, typename, int> class Evaluator>
class PoissonSeries {
 public:
  using AperiodicPolynomial =
      PolynomialInMonomialBasis<Value, Instant, aperiodic_degree_, Evaluator>;
  using PeriodicPolynomial =
      PolynomialInMonomialBasis<Value, Instant, periodic_degree_, Evaluator>;

  // TODO(phl): Use designated initializers for this struct once this project
  // can be compiled using c++latest.
  // TODO(phl): If we wanted to have Poisson series returning affine values,
  // these polynomials should be changed to return Difference<Value>.
  struct Polynomials {
    PeriodicPolynomial sin;
    PeriodicPolynomial cos;
  };

  using PolynomialsByAngularFrequency =
      std::vector<std::pair<AngularFrequency, Polynomials>>;

  // The |periodic| vector may contain frequencies that are negative or zero, as
  // well as repeated frequencies.  It is not expected to be ordered.
  PoissonSeries(AperiodicPolynomial const& aperiodic,
                PolynomialsByAngularFrequency const& periodic);

  // A Poisson series may be explicitly converted to a higher degree (possibly
  // with a different evaluator).
  template<int higher_aperiodic_degree, int higher_periodic_degree,
           template<typename, typename, int> class HigherEvaluator>
  explicit operator PoissonSeries<Value,
                                  higher_aperiodic_degree,
                                  higher_periodic_degree,
                                  HigherEvaluator>() const;

  Instant const& origin() const;

  Value operator()(Instant const& t) const;

  // Returns a copy of this series adjusted to the given origin.
  PoissonSeries AtOrigin(Instant const& origin) const;

  // The constant term of the result is zero.
  PoissonSeries<quantities::Primitive<Value, Time>,
                aperiodic_degree_ + 1, periodic_degree_ + 1,
                Evaluator>
  Primitive() const;

  quantities::Primitive<Value, Time> Integrate(Instant const& t1,
                                               Instant const& t2) const;

  template<int aperiodic_wdegree, int periodic_wdegree>
  typename Hilbert<Value>::NormType Norm(
      PoissonSeries<double,
                    aperiodic_wdegree, periodic_wdegree,
                    Evaluator> const& weight,
                    Instant const& t_min,
                    Instant const& t_max) const;

  template<int aperiodic_rdegree, int periodic_rdegree>
  PoissonSeries& operator+=(PoissonSeries<Value,
                                          aperiodic_rdegree, periodic_rdegree,
                                          Evaluator> const& right);
  template<int aperiodic_rdegree, int periodic_rdegree>
  PoissonSeries& operator-=(PoissonSeries<Value,
                                          aperiodic_rdegree, periodic_rdegree,
                                          Evaluator> const& right);

  void WriteToMessage(not_null<serialization::PoissonSeries*> message) const;
  static PoissonSeries ReadFromMessage(
      serialization::PoissonSeries const& message);

 private:
  // Similar to the public constructor, but passing by copy allows moves, which
  // is useful for internal algorithms.
  struct PrivateConstructor {};
  PoissonSeries(PrivateConstructor,
                AperiodicPolynomial aperiodic,
                PolynomialsByAngularFrequency periodic);

  // Similar to the previous constructor, except that the |periodic| vector is
  // used verbatim, without sorting or normalization, which is useful for
  // internal algorithms which produce positive, ordered frequencies.
  struct TrustedPrivateConstructor {};
  PoissonSeries(TrustedPrivateConstructor,
                AperiodicPolynomial aperiodic,
                PolynomialsByAngularFrequency periodic);

  // Splits this series into two copies, with frequencies lower and higher than
  // ω_cutoff, respectively.
  struct SplitPoissonSeries {
    PoissonSeries slow;
    PoissonSeries fast;
  };
  SplitPoissonSeries Split(AngularFrequency const& ω_cutoff) const;

  Instant origin_;  // Common to all polynomials.
  AperiodicPolynomial aperiodic_;
  // The frequencies in this vector are positive, distinct and in increasing
  // order.
  PolynomialsByAngularFrequency periodic_;

  template<typename V, int ar, int pr,
           template<typename, typename, int> class E>
  PoissonSeries<V, ar, pr, E>
  friend operator-(PoissonSeries<V, ar, pr, E> const& right);
  template<typename V, int al, int pl, int ar, int pr,
           template<typename, typename, int> class E>
  PoissonSeries<V, PRINCIPIA_MAX(al, ar), PRINCIPIA_MAX(pl, pr), E>
  friend operator+(PoissonSeries<V, al, pl, E> const& left,
                   PoissonSeries<V, ar, pr, E> const& right);
  template<typename V, int al, int pl, int ar, int pr,
           template<typename, typename, int> class E>
  PoissonSeries<V, PRINCIPIA_MAX(al, ar), PRINCIPIA_MAX(pl, pr), E>
  friend operator-(PoissonSeries<V, al, pl, E> const& left,
                   PoissonSeries<V, ar, pr, E> const& right);
  template<typename S, typename V, int ar, int pr,
           template<typename, typename, int> class E>
  PoissonSeries<Product<S, V>, ar, pr, E>
  friend operator*(S const& left,
                   PoissonSeries<V, ar, pr, E> const& right);
  template<typename S, typename V, int al, int pl,
           template<typename, typename, int> class E>
  PoissonSeries<Product<V, S>, al, pl, E>
  friend operator*(PoissonSeries<V, al, pl, E> const& left,
                   S const& right);
  template<typename S, typename V, int al, int pl,
           template<typename, typename, int> class E>
  PoissonSeries<Quotient<V, S>, al, pl, E>
  friend operator/(PoissonSeries<V, al, pl, E> const& left,
                   S const& right);
  template<typename L, typename R, int al, int pl, int ar, int pr,
           template<typename, typename, int> class E,
           typename P>
  auto friend Multiply(PoissonSeries<L, al, pl, E> const& left,
                       PoissonSeries<R, ar, pr, E> const& right,
                       P const& product);
  template<typename V, int ad, int pd,
           template<typename, typename, int> class E>
  friend std::ostream& operator<<(std::ostream& out,
                                  PoissonSeries<V, ad, pd, E> const& series);
  template<typename L, typename R,
         int al, int pl, int ar, int pr, int aw, int pw,
         template<typename, typename, int> class E>
  friend typename Hilbert<L, R>::InnerProductType InnerProduct(
      PoissonSeries<L, al, pl, E> const& left,
      PoissonSeries<R, ar, pr, E> const& right,
      PoissonSeries<double, aw, pw, E> const& weight,
      Instant const& t_min,
      Instant const& t_max);
  template<typename V, int ad, int pd,
           template<typename, typename, int> class E,
           typename O>
  friend std::string mathematica::internal_mathematica::ToMathematicaExpression(
      PoissonSeries<V, ad, pd, E> const& polynomial,
      O express_in);
  template<typename V, int ad, int pd,
           template<typename, typename, int> class E>
  friend class PoissonSeries;
};

// Vector space of Poisson series.

template<typename Value,
         int aperiodic_rdegree, int periodic_rdegree,
         template<typename, typename, int> class Evaluator>
PoissonSeries<Value, aperiodic_rdegree, periodic_rdegree, Evaluator>
operator+(PoissonSeries<Value,
                        aperiodic_rdegree, periodic_rdegree,
                        Evaluator> const& right);

template<typename Value,
         int aperiodic_rdegree, int periodic_rdegree,
         template<typename, typename, int> class Evaluator>
PoissonSeries<Value, aperiodic_rdegree, periodic_rdegree, Evaluator>
operator-(PoissonSeries<Value,
                        aperiodic_rdegree, periodic_rdegree,
                        Evaluator> const& right);

template<typename Value,
         int aperiodic_ldegree, int periodic_ldegree,
         int aperiodic_rdegree, int periodic_rdegree,
         template<typename, typename, int> class Evaluator>
PoissonSeries<Value,
              PRINCIPIA_MAX(aperiodic_ldegree, aperiodic_rdegree),
              PRINCIPIA_MAX(periodic_ldegree, periodic_rdegree),
              Evaluator>
operator+(PoissonSeries<Value,
                        aperiodic_ldegree, periodic_ldegree,
                        Evaluator> const& left,
          PoissonSeries<Value,
                        aperiodic_rdegree, periodic_rdegree,
                        Evaluator> const& right);

template<typename Value,
         int aperiodic_ldegree, int periodic_ldegree,
         int aperiodic_rdegree, int periodic_rdegree,
         template<typename, typename, int> class Evaluator>
PoissonSeries<Value,
              PRINCIPIA_MAX(aperiodic_ldegree, aperiodic_rdegree),
              PRINCIPIA_MAX(periodic_ldegree, periodic_rdegree),
              Evaluator>
operator-(PoissonSeries<Value,
                        aperiodic_ldegree, periodic_ldegree,
                        Evaluator> const& left,
          PoissonSeries<Value,
                        aperiodic_rdegree, periodic_rdegree,
                        Evaluator> const& right);

template<typename Scalar,
         typename Value,
         int aperiodic_rdegree, int periodic_rdegree,
         template<typename, typename, int> class Evaluator>
PoissonSeries<Product<Scalar, Value>,
              aperiodic_rdegree, periodic_rdegree,
              Evaluator>
operator*(Scalar const& left,
          PoissonSeries<Value,
                        aperiodic_rdegree, periodic_rdegree,
                        Evaluator> const& right);

template<typename Scalar,
         typename Value,
         int aperiodic_ldegree, int periodic_ldegree,
         template<typename, typename, int> class Evaluator>
PoissonSeries<Product<Value, Scalar>,
              aperiodic_ldegree, periodic_ldegree,
              Evaluator>
operator*(PoissonSeries<Value,
                        aperiodic_ldegree, periodic_ldegree,
                        Evaluator> const& left,
          Scalar const& right);

template<typename Scalar,
         typename Value,
         int aperiodic_ldegree, int periodic_ldegree,
         template<typename, typename, int> class Evaluator>
PoissonSeries<Quotient<Value, Scalar>,
              aperiodic_ldegree, periodic_ldegree,
              Evaluator>
operator/(PoissonSeries<Value,
                        aperiodic_ldegree, periodic_ldegree,
                        Evaluator> const& left,
          Scalar const& right);

// Algebra of Poisson series.

template<typename LValue, typename RValue,
         int aperiodic_ldegree, int periodic_ldegree,
         int aperiodic_rdegree, int periodic_rdegree,
         template<typename, typename, int> class Evaluator>
PoissonSeries<Product<LValue, RValue>,
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
          PoissonSeries<RValue,
                        aperiodic_rdegree, periodic_rdegree,
                        Evaluator> const& right);

// Returns a scalar-valued Poisson series obtained by pointwise inner product of
// two vector-valued series.
template<typename LValue, typename RValue,
         int aperiodic_ldegree, int periodic_ldegree,
         int aperiodic_rdegree, int periodic_rdegree,
         template<typename, typename, int> class Evaluator>
PoissonSeries<typename Hilbert<LValue, RValue>::InnerProductType,
              PRINCIPIA_MAX4(aperiodic_ldegree + aperiodic_rdegree,
                             aperiodic_ldegree + periodic_rdegree,
                             periodic_ldegree + aperiodic_rdegree,
                             periodic_ldegree + periodic_rdegree),
              PRINCIPIA_MAX3(aperiodic_ldegree + periodic_rdegree,
                             periodic_ldegree + aperiodic_rdegree,
                             periodic_ldegree + periodic_rdegree),
              Evaluator>
PointwiseInnerProduct(PoissonSeries<LValue,
                                    aperiodic_ldegree, periodic_ldegree,
                                    Evaluator> const& left,
                      PoissonSeries<RValue,
                                    aperiodic_rdegree, periodic_rdegree,
                                    Evaluator> const& right);

// Output.

template<typename Value,
         int aperiodic_degree, int periodic_degree,
         template<typename, typename, int> class Evaluator>
std::ostream& operator<<(
    std::ostream& out,
    PoissonSeries<Value,
                  aperiodic_degree, periodic_degree,
                  Evaluator> const& series);

// Inner product space of Poisson series.

// Technically the weight function must be nonnegative for this to be an inner
// product.  Not sure how this works with the flat-top windows, which can be
// negative.  Note that the result is normalized by dividing by (t_max - t_min).
template<typename LValue, typename RValue,
         int aperiodic_ldegree, int periodic_ldegree,
         int aperiodic_rdegree, int periodic_rdegree,
         int aperiodic_wdegree, int periodic_wdegree,
         template<typename, typename, int> class Evaluator>
typename Hilbert<LValue, RValue>::InnerProductType
InnerProduct(PoissonSeries<LValue,
                           aperiodic_ldegree, periodic_ldegree,
                           Evaluator> const& left,
             PoissonSeries<RValue,
                           aperiodic_rdegree, periodic_rdegree,
                           Evaluator> const& right,
             PoissonSeries<double,
                           aperiodic_wdegree, periodic_wdegree,
                           Evaluator> const& weight,
             Instant const& t_min,
             Instant const& t_max);

}  // namespace internal_poisson_series

using internal_poisson_series::InnerProduct;
using internal_poisson_series::PoissonSeries;

}  // namespace numerics
}  // namespace principia

#include "numerics/poisson_series_body.hpp"
