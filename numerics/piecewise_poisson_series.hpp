
#pragma once

#include <algorithm>
#include <optional>
#include <string>
#include <vector>

#include "base/macros.hpp"
#include "base/not_null.hpp"
#include "geometry/complexification.hpp"
#include "geometry/hilbert.hpp"
#include "geometry/interval.hpp"
#include "geometry/named_quantities.hpp"
#include "numerics/poisson_series.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"
#include "serialization/numerics.pb.h"

#if PRINCIPIA_COMPILER_MSVC_HAS_CXX20
#define PRINCIPIA_MAX(l, r) ((l) > (r) ? (l) : (r))
#define PRINCIPIA_MAX3(x1, x2, x3) \
  PRINCIPIA_MAX((x1), PRINCIPIA_MAX((x2), (x3)))
#define PRINCIPIA_MAX4(x1, x2, x3, x4) \
  PRINCIPIA_MAX((x1), PRINCIPIA_MAX((x2), PRINCIPIA_MAX((x3), (x4))))
#else
#define PRINCIPIA_MAX(l, r) std::max((l), (r))
#define PRINCIPIA_MAX3(x1, x2, x3) std::max({(x1), (x2), (x3)})
#define PRINCIPIA_MAX4(x1, x2, x3) std::max({(x1), (x2), (x3), (x4)})
#endif

namespace principia {
namespace numerics {
FORWARD_DECLARE_FROM(piecewise_poisson_series,
                     TEMPLATE(typename Value,
                              int aperiodic_degree, int periodic_degree,
                              template<typename, typename, int> class Evaluator)
                              class,
                     PiecewisePoissonSeries);
}  // namespace numerics

namespace mathematica {
FORWARD_DECLARE_FUNCTION_FROM(
    mathematica,
    TEMPLATE(typename Value,
             int aperiodic_degree, int periodic_degree,
             template<typename, typename, int> class Evaluator,
             typename OptionalExpressIn) std::string,
    ToMathematicaExpression,
    (numerics::PiecewisePoissonSeries<Value,
                                      aperiodic_degree, periodic_degree,
                                      Evaluator> const& series,
     OptionalExpressIn express_in));
}  // namespace mathematica

namespace numerics {
namespace internal_piecewise_poisson_series {

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

// The trigonometric functions are by default assumed to look like a polynomial
// of this degree over an interval of a piecewise series.
constexpr int estimated_trigonometric_degree = 1;

// A function defined by Poisson series piecewise.  Each of the Poisson series
// making up the function applies over the semi-open interval
// [internal.min, interval.max[.  It's not required that the function be
// continuous.
template<typename Value,
         int aperiodic_degree_, int periodic_degree_,
         template<typename, typename, int> class Evaluator>
class PiecewisePoissonSeries {
 public:
  using Series = PoissonSeries<Value,
                               aperiodic_degree_, periodic_degree_,
                               Evaluator>;
  using Spectrum = std::function<Complexification<Primitive<Value, Time>>(
      AngularFrequency const&)>;

  PiecewisePoissonSeries(Interval<Instant> const& interval,
                         Series const& series);

  // The intervals for successive calls to Append must be consecutive.  For the
  // first call, the interval must be consecutive with the one passed at
  // construction.
  void Append(Interval<Instant> const& interval,
              Series const& series);

  Instant t_min() const;
  Instant t_max() const;

  // t must be in the interval [t_min, t_max].
  Value operator()(Instant const& t) const;

  // Returns the Fourier transform of this piecewise Poisson series.
  // The function is taken to be 0 outside [t_min, t_max].
  // The convention used is ∫ f(t) exp(-iωt) dt, corresponding to Mathematica’s
  // FourierParameters -> {1, -1} for FourierTransform (the “pure mathematics;
  // systems engineering”).
  // When evaluated at a given frequency, the Fourier transform is computed by
  // Gauss-Legendre quadrature on each subinterval, where the number of points
  // is chosen assuming that the periods of periodic terms are all large
  // compared to the subintervals.
  // If apodization is desired, |*this| should be multiplied by an apodization
  // function, and |FourierTransform| should be called on the product.
  // |*this| must outlive the resulting function.
  Spectrum FourierTransform() const;

  template<int aperiodic_rdegree, int periodic_rdegree>
  PiecewisePoissonSeries& operator+=(
      PoissonSeries<Value,
                    aperiodic_rdegree, periodic_rdegree,
                    Evaluator> const& right);
  template<int aperiodic_rdegree, int periodic_rdegree>
  PiecewisePoissonSeries& operator-=(
      PoissonSeries<Value,
                    aperiodic_rdegree, periodic_rdegree,
                    Evaluator> const& right);

  void WriteToMessage(
      not_null<serialization::PiecewisePoissonSeries*> message) const;
  static PiecewisePoissonSeries ReadFromMessage(
      serialization::PiecewisePoissonSeries const& message);

 private:
  PiecewisePoissonSeries(std::vector<Instant> const& bounds,
                         std::vector<Series> const& series,
                         std::optional<Series> const& addend);

  Value EvaluateAddend(Instant const& t) const;

  std::vector<Instant> bounds_;
  std::vector<Series> series_;
  std::optional<Series> addend_;

  template<typename V, int ar, int pr,
           template<typename, typename, int> class E>
  PiecewisePoissonSeries<V, ar, pr, E>
  friend operator-(PiecewisePoissonSeries<V, ar, pr, E> const& right);
  template<typename S, typename V, int ar, int pr,
           template<typename, typename, int> class E>
  PiecewisePoissonSeries<Product<S, V>, ar, pr, E>
  friend operator*(S const& left,
                   PiecewisePoissonSeries<V, ar, pr, E> const& right);
  template<typename S, typename V, int al, int pl,
           template<typename, typename, int> class E>
  PiecewisePoissonSeries<Product<V, S>, al, pl, E>
  friend operator*(PiecewisePoissonSeries<V, al, pl, E> const& left,
                   S const& right);
  template<typename S, typename V, int al, int pl,
           template<typename, typename, int> class E>
  PiecewisePoissonSeries<Quotient<V, S>, al, pl, E>
  friend operator/(PiecewisePoissonSeries<V, al, pl, E> const& left,
                   S const& right);
  template<typename V, int al, int pl, int ar, int pr,
           template<typename, typename, int> class E>
  PiecewisePoissonSeries<V, PRINCIPIA_MAX(al, ar), PRINCIPIA_MAX(pl, pr), E>
  friend operator+(PoissonSeries<V, al, pl, E> const& left,
                   PiecewisePoissonSeries<V, ar, pr, E> const& right);
  template<typename V, int al, int pl, int ar, int pr,
           template<typename, typename, int> class E>
  PiecewisePoissonSeries<V, PRINCIPIA_MAX(al, ar), PRINCIPIA_MAX(pl, pr), E>
  friend operator+(PiecewisePoissonSeries<V, al, pl, E> const& left,
                   PoissonSeries<V, ar, pr, E> const& right);
  template<typename V, int al, int pl, int ar, int pr,
           template<typename, typename, int> class E>
  PiecewisePoissonSeries<V, PRINCIPIA_MAX(al, ar), PRINCIPIA_MAX(pl, pr), E>
  friend operator-(PoissonSeries<V, al, pl, E> const& left,
                   PiecewisePoissonSeries<V, ar, pr, E> const& right);
  template<typename V, int al, int pl, int ar, int pr,
           template<typename, typename, int> class E>
  PiecewisePoissonSeries<V, PRINCIPIA_MAX(al, ar), PRINCIPIA_MAX(pl, pr), E>
  friend operator-(PiecewisePoissonSeries<V, al, pl, E> const& left,
                   PoissonSeries<V, ar, pr, E> const& right);
  template<typename L, typename R, int al, int pl, int ar, int pr,
           template<typename, typename, int> class E>
  PiecewisePoissonSeries<Product<L, R>,
                         PRINCIPIA_MAX4(al + ar, al + pr, pl + ar, pl + pr),
                         PRINCIPIA_MAX3(al + pr, pl + ar, pl + pr),
                         E>
  friend operator*(PoissonSeries<L, al, pl, E> const& left,
                   PiecewisePoissonSeries<R, ar, pr, E> const& right);
  template<typename L, typename R, int al, int pl, int ar, int pr,
           template<typename, typename, int> class E>
  PiecewisePoissonSeries<Product<L, R>,
                         PRINCIPIA_MAX4(al + ar, al + pr, pl + ar, pl + pr),
                         PRINCIPIA_MAX3(al + pr, pl + ar, pl + pr),
                         E>
  friend operator*(PiecewisePoissonSeries<L, al, pl, E> const& left,
                   PoissonSeries<R, ar, pr, E> const& right);
  template<typename L, typename R,
           int al, int pl, int ar, int pr, int aw, int pw,
           template<typename, typename, int> class E, int p>
  typename Hilbert<L, R>::InnerProductType
  friend InnerProduct(PoissonSeries<L, al, pl, E> const& left,
                      PiecewisePoissonSeries<R, ar, pr, E> const& right,
                      PoissonSeries<double, aw, pw, E> const& weight,
                      Instant const& t_min,
                      Instant const& t_max);
  template<typename L, typename R,
           int al, int pl, int ar, int pr, int aw, int pw,
           template<typename, typename, int> class E, int p>
  typename Hilbert<L, R>::InnerProductType
  friend InnerProduct(PiecewisePoissonSeries<L, al, pl, E> const& left,
                      PoissonSeries<R, ar, pr, E> const& right,
                      PoissonSeries<double, aw, pw, E> const& weight,
                      Instant const& t_min,
                      Instant const& t_max);
  template<typename V, int ad, int pd,
           template<typename, typename, int> class E,
           typename O>
  friend std::string mathematica::internal_mathematica::ToMathematicaExpression(
      PiecewisePoissonSeries<V, ad, pd, E> const& polynomial,
      O express_in);
};

// Some of the vector space operations for piecewise Poisson series.  Note that
// while piecewise Poisson series have an algebra structure, we do not need it.

template<typename Value,
         int aperiodic_rdegree, int periodic_rdegree,
         template<typename, typename, int> class Evaluator>
PiecewisePoissonSeries<Value, aperiodic_rdegree, periodic_rdegree, Evaluator>
operator+(PiecewisePoissonSeries<Value,
                                 aperiodic_rdegree, periodic_rdegree,
                                 Evaluator> const& right);

template<typename Value,
         int aperiodic_rdegree, int periodic_rdegree,
         template<typename, typename, int> class Evaluator>
PiecewisePoissonSeries<Value, aperiodic_rdegree, periodic_rdegree, Evaluator>
operator-(PiecewisePoissonSeries<Value,
                                 aperiodic_rdegree, periodic_rdegree,
                                 Evaluator> const& right);

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
                                 Evaluator> const& right);

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
          Scalar const& right);

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
          Scalar const& right);

// Action of Poisson series on piecewise Poisson series.

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
                                 Evaluator> const& right);

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
                        Evaluator> const& right);

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
                                 Evaluator> const& right);

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
                        Evaluator> const& right);

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
                                 Evaluator> const& right);

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
                        Evaluator> const& right);

template<
    typename LValue, typename RValue,
    int aperiodic_ldegree, int periodic_ldegree,
    int aperiodic_rdegree, int periodic_rdegree,
    int aperiodic_wdegree, int periodic_wdegree,
    template<typename, typename, int> class Evaluator,
    int points = (std::max(aperiodic_ldegree,
                           periodic_ldegree + estimated_trigonometric_degree) +
                  std::max(aperiodic_ldegree,
                           periodic_ldegree + estimated_trigonometric_degree) +
                  std::max(aperiodic_ldegree,
                           periodic_ldegree + estimated_trigonometric_degree)) /
                 2>
typename Hilbert<LValue, RValue>::InnerProductType InnerProduct(
    PoissonSeries<LValue,
                  aperiodic_ldegree, periodic_ldegree,
                  Evaluator> const& left,
    PiecewisePoissonSeries<RValue,
                           aperiodic_rdegree, periodic_rdegree,
                           Evaluator> const& right,
    PoissonSeries<double,
                  aperiodic_wdegree, periodic_wdegree,
                  Evaluator> const& weight);

template<
    typename LValue, typename RValue,
    int aperiodic_ldegree, int periodic_ldegree,
    int aperiodic_rdegree, int periodic_rdegree,
    int aperiodic_wdegree, int periodic_wdegree,
    template<typename, typename, int> class Evaluator,
    int points = (std::max(aperiodic_ldegree,
                           periodic_ldegree + estimated_trigonometric_degree) +
                  std::max(aperiodic_ldegree,
                           periodic_ldegree + estimated_trigonometric_degree) +
                  std::max(aperiodic_ldegree,
                           periodic_ldegree + estimated_trigonometric_degree)) /
                 2>
typename Hilbert<LValue, RValue>::InnerProductType InnerProduct(
    PoissonSeries<LValue,
                  aperiodic_ldegree, periodic_ldegree,
                  Evaluator> const& left,
    PiecewisePoissonSeries<RValue,
                           aperiodic_rdegree, periodic_rdegree,
                           Evaluator> const& right,
    PoissonSeries<double,
                  aperiodic_wdegree, periodic_wdegree,
                  Evaluator> const& weight,
    Instant const& t_min,
    Instant const& t_max);

template<
    typename LValue, typename RValue,
    int aperiodic_ldegree, int periodic_ldegree,
    int aperiodic_rdegree, int periodic_rdegree,
    int aperiodic_wdegree, int periodic_wdegree,
    template<typename, typename, int> class Evaluator,
    int points = (std::max(aperiodic_ldegree,
                           periodic_ldegree + estimated_trigonometric_degree) +
                  std::max(aperiodic_ldegree,
                           periodic_ldegree + estimated_trigonometric_degree) +
                  std::max(aperiodic_ldegree,
                           periodic_ldegree + estimated_trigonometric_degree)) /
                 2>
typename Hilbert<LValue, RValue>::InnerProductType InnerProduct(
    PiecewisePoissonSeries<LValue,
                           aperiodic_ldegree, periodic_ldegree,
                           Evaluator> const& left,
    PoissonSeries<RValue,
                  aperiodic_rdegree, periodic_rdegree,
                  Evaluator> const& right,
    PoissonSeries<double,
                  aperiodic_wdegree, periodic_wdegree,
                  Evaluator> const& weight);

template<
    typename LValue, typename RValue,
    int aperiodic_ldegree, int periodic_ldegree,
    int aperiodic_rdegree, int periodic_rdegree,
    int aperiodic_wdegree, int periodic_wdegree,
    template<typename, typename, int> class Evaluator,
    int points = (std::max(aperiodic_ldegree,
                           periodic_ldegree + estimated_trigonometric_degree) +
                  std::max(aperiodic_ldegree,
                           periodic_ldegree + estimated_trigonometric_degree) +
                  std::max(aperiodic_ldegree,
                           periodic_ldegree + estimated_trigonometric_degree)) /
                 2>
typename Hilbert<LValue, RValue>::InnerProductType InnerProduct(
    PiecewisePoissonSeries<LValue,
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

}  // namespace internal_piecewise_poisson_series

using internal_piecewise_poisson_series::PiecewisePoissonSeries;

}  // namespace numerics
}  // namespace principia

#include "numerics/piecewise_poisson_series_body.hpp"

#undef PRINCIPIA_MAX
#undef PRINCIPIA_MAX3
#undef PRINCIPIA_MAX4
