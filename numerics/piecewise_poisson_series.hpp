#pragma once

#include <algorithm>
#include <memory>
#include <optional>
#include <string>
#include <vector>

#include "base/macros.hpp"  // 🧙 For forward declarations.
#include "base/not_null.hpp"
#include "geometry/complexification.hpp"
#include "geometry/hilbert.hpp"
#include "geometry/instant.hpp"
#include "geometry/interval.hpp"
#include "numerics/poisson_series.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"
#include "serialization/numerics.pb.h"

namespace principia {
namespace numerics {
FORWARD_DECLARE(
    TEMPLATE(typename Value,
             int aperiodic_degree, int periodic_degree) class,
    PiecewisePoissonSeries,
    FROM(piecewise_poisson_series));
}  // namespace numerics

namespace mathematica {
FORWARD_DECLARE_FUNCTION(
    TEMPLATE(typename Value,
             int aperiodic_degree, int periodic_degree,
             typename OptionalExpressIn) std::string,
    ToMathematicaBody,
    (numerics::_piecewise_poisson_series::PiecewisePoissonSeries<
         Value,
         aperiodic_degree, periodic_degree> const& series,
     OptionalExpressIn express_in),
     FROM(mathematica));
}  // namespace mathematica

namespace numerics {
namespace _piecewise_poisson_series {
namespace internal {

using namespace principia::base::_not_null;
using namespace principia::geometry::_complexification;
using namespace principia::geometry::_hilbert;
using namespace principia::geometry::_instant;
using namespace principia::geometry::_interval;
using namespace principia::numerics::_poisson_series;
using namespace principia::quantities::_named_quantities;
using namespace principia::quantities::_quantities;

// The trigonometric functions are by default assumed to look like a polynomial
// of this degree over an interval of a piecewise series.
constexpr int estimated_trigonometric_degree = 1;

// A function defined by Poisson series piecewise.  Each of the Poisson series
// making up the function applies over the semi-open interval
// [internal.min, interval.max[.  It's not required that the function be
// continuous.
template<typename Value,
         int aperiodic_degree_, int periodic_degree_>
class PiecewisePoissonSeries {
 public:
  using Series = PoissonSeries<Value,
                               aperiodic_degree_, periodic_degree_>;
  using Spectrum = std::function<Complexification<Primitive<Value, Time>>(
      AngularFrequency const&)>;

  PiecewisePoissonSeries(Interval<Instant> const& interval,
                         Series const& series);

  template<int aperiodic_rdegree, int periodic_rdegree>
  PiecewisePoissonSeries& operator+=(
      PoissonSeries<Value,
                    aperiodic_rdegree, periodic_rdegree> const& right);
  template<int aperiodic_rdegree, int periodic_rdegree>
  PiecewisePoissonSeries& operator-=(
      PoissonSeries<Value,
                    aperiodic_rdegree, periodic_rdegree> const& right);

  // The intervals for successive calls to Append must be consecutive.  For the
  // first call, the interval must be consecutive with the one passed at
  // construction.
  void Append(Interval<Instant> const& interval,
              Series const& series);

  Instant t_min() const;
  Instant t_max() const;
  AngularFrequency max_ω() const;

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
  // If apodization is desired, `*this` should be multiplied by an apodization
  // function, and `FourierTransform` should be called on the product.
  // `*this` must outlive the resulting function.
  Spectrum FourierTransform() const;

  template<int aperiodic_wdegree, int periodic_wdegree>
  typename Hilbert<Value>::NormType Norm(
      PoissonSeries<double,
                    aperiodic_wdegree, periodic_wdegree> const& weight,
      Instant const& t_min,
      Instant const& t_max) const;

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

  // A cache of the evaluation of this function at the Gauss-Legendre points.
  // Caching is possible because the Gauss-Legendre integration picks its points
  // deterministically.  Note that caching keeps working if new series get
  // appended to this function: the cache will just be "too short" and new
  // entries can be appended as needed.  We use a shared_ptr<> to make sure that
  // copies remain cheap.
  mutable std::shared_ptr<std::vector<Value>> gauss_legendre_cache_;

  template<typename V, int ar, int pr>
  PiecewisePoissonSeries<V, ar, pr>
  friend operator-(PiecewisePoissonSeries<V, ar, pr> const& right);
  template<typename S, typename V, int ar, int pr>
  PiecewisePoissonSeries<Product<S, V>, ar, pr>
  friend operator*(S const& left,
                   PiecewisePoissonSeries<V, ar, pr> const& right);
  template<typename S, typename V, int al, int pl>
  PiecewisePoissonSeries<Product<V, S>, al, pl>
  friend operator*(PiecewisePoissonSeries<V, al, pl> const& left,
                   S const& right);
  template<typename S, typename V, int al, int pl>
  PiecewisePoissonSeries<Quotient<V, S>, al, pl>
  friend operator/(PiecewisePoissonSeries<V, al, pl> const& left,
                   S const& right);
  template<typename V, int al, int pl, int ar, int pr>
  PiecewisePoissonSeries<V, std::max(al, ar), std::max(pl, pr)>
  friend operator+(PoissonSeries<V, al, pl> const& left,
                   PiecewisePoissonSeries<V, ar, pr> const& right);
  template<typename V, int al, int pl, int ar, int pr>
  PiecewisePoissonSeries<V, std::max(al, ar), std::max(pl, pr)>
  friend operator+(PiecewisePoissonSeries<V, al, pl> const& left,
                   PoissonSeries<V, ar, pr> const& right);
  template<typename V, int al, int pl, int ar, int pr>
  PiecewisePoissonSeries<V, std::max(al, ar), std::max(pl, pr)>
  friend operator-(PoissonSeries<V, al, pl> const& left,
                   PiecewisePoissonSeries<V, ar, pr> const& right);
  template<typename V, int al, int pl, int ar, int pr>
  PiecewisePoissonSeries<V, std::max(al, ar), std::max(pl, pr)>
  friend operator-(PiecewisePoissonSeries<V, al, pl> const& left,
                   PoissonSeries<V, ar, pr> const& right);
  template<typename L, typename R, int al, int pl, int ar, int pr>
  PiecewisePoissonSeries<Product<L, R>,
                         std::max({al + ar, al + pr, pl + ar, pl + pr}),
                         std::max({al + pr, pl + ar, pl + pr})>
  friend operator*(PoissonSeries<L, al, pl> const& left,
                   PiecewisePoissonSeries<R, ar, pr> const& right);
  template<typename L, typename R, int al, int pl, int ar, int pr>
  PiecewisePoissonSeries<Product<L, R>,
                         std::max({al + ar, al + pr, pl + ar, pl + pr}),
                         std::max({al + pr, pl + ar, pl + pr})>
  friend operator*(PiecewisePoissonSeries<L, al, pl> const& left,
                   PoissonSeries<R, ar, pr> const& right);
  template<typename L, typename R,
           int al, int pl, int ar, int pr, int aw, int pw>
  typename Hilbert<L, R>::InnerProductType
  friend InnerProduct(PoissonSeries<L, al, pl> const& left,
                      PiecewisePoissonSeries<R, ar, pr> const& right,
                      PoissonSeries<double, aw, pw> const& weight,
                      Instant const& t_min,
                      Instant const& t_max,
                      std::optional<int> max_points);
  template<typename L, typename R,
           int al, int pl, int ar, int pr, int aw, int pw>
  typename Hilbert<L, R>::InnerProductType
  friend InnerProduct(PiecewisePoissonSeries<L, al, pl> const& left,
                      PoissonSeries<R, ar, pr> const& right,
                      PoissonSeries<double, aw, pw> const& weight,
                      Instant const& t_min,
                      Instant const& t_max,
                      std::optional<int> max_points);
  template<typename V, int ad, int pd, typename O>
  friend std::string mathematica::_mathematica::internal::ToMathematicaBody(
      PiecewisePoissonSeries<V, ad, pd> const& polynomial,
      O express_in);
};

// Some of the vector space operations for piecewise Poisson series.  Note that
// while piecewise Poisson series have an algebra structure, we do not need it.

template<typename Value,
         int aperiodic_rdegree, int periodic_rdegree>
PiecewisePoissonSeries<Value, aperiodic_rdegree, periodic_rdegree>
operator+(PiecewisePoissonSeries<
              Value, aperiodic_rdegree, periodic_rdegree> const& right);

template<typename Value,
         int aperiodic_rdegree, int periodic_rdegree>
PiecewisePoissonSeries<Value, aperiodic_rdegree, periodic_rdegree>
operator-(PiecewisePoissonSeries<
               Value, aperiodic_rdegree, periodic_rdegree> const& right);

template<typename Scalar,
         typename Value,
         int aperiodic_rdegree, int periodic_rdegree>
PiecewisePoissonSeries<Product<Scalar, Value>,
                       aperiodic_rdegree, periodic_rdegree>
operator*(Scalar const& left,
          PiecewisePoissonSeries<
              Value, aperiodic_rdegree, periodic_rdegree> const& right);

template<typename Scalar,
         typename Value,
         int aperiodic_ldegree, int periodic_ldegree>
PiecewisePoissonSeries<Product<Value, Scalar>,
                       aperiodic_ldegree, periodic_ldegree>
operator*(PiecewisePoissonSeries<
              Value, aperiodic_ldegree, periodic_ldegree> const& left,
          Scalar const& right);

template<typename Scalar,
         typename Value,
         int aperiodic_ldegree, int periodic_ldegree>
PiecewisePoissonSeries<Quotient<Value, Scalar>,
                       aperiodic_ldegree, periodic_ldegree>
operator/(PiecewisePoissonSeries<
              Value, aperiodic_ldegree, periodic_ldegree> const& left,
          Scalar const& right);

// Action of Poisson series on piecewise Poisson series.

template<typename Value,
         int aperiodic_ldegree, int periodic_ldegree,
         int aperiodic_rdegree, int periodic_rdegree>
PiecewisePoissonSeries<Value,
                       std::max(aperiodic_ldegree, aperiodic_rdegree),
                       std::max(periodic_ldegree, periodic_rdegree)>
operator+(PoissonSeries<Value,
                        aperiodic_ldegree, periodic_ldegree> const& left,
          PiecewisePoissonSeries<
               Value, aperiodic_rdegree, periodic_rdegree> const& right);

template<typename Value,
         int aperiodic_ldegree, int periodic_ldegree,
         int aperiodic_rdegree, int periodic_rdegree>
PiecewisePoissonSeries<Value,
                       std::max(aperiodic_ldegree, aperiodic_rdegree),
                       std::max(periodic_ldegree, periodic_rdegree)>
operator+(PiecewisePoissonSeries<
              Value, aperiodic_ldegree, periodic_ldegree> const& left,
          PoissonSeries<Value,
                        aperiodic_rdegree, periodic_rdegree> const& right);

template<typename Value,
         int aperiodic_ldegree, int periodic_ldegree,
         int aperiodic_rdegree, int periodic_rdegree>
PiecewisePoissonSeries<Value,
                       std::max(aperiodic_ldegree, aperiodic_rdegree),
                       std::max(periodic_ldegree, periodic_rdegree)>
operator-(PoissonSeries<Value,
                        aperiodic_ldegree, periodic_ldegree> const& left,
          PiecewisePoissonSeries<
              Value, aperiodic_rdegree, periodic_rdegree> const& right);

template<typename Value,
         int aperiodic_ldegree, int periodic_ldegree,
         int aperiodic_rdegree, int periodic_rdegree>
PiecewisePoissonSeries<Value,
                       std::max(aperiodic_ldegree, aperiodic_rdegree),
                       std::max(periodic_ldegree, periodic_rdegree)>
operator-(PiecewisePoissonSeries<
              Value, aperiodic_ldegree, periodic_ldegree> const& left,
          PoissonSeries<Value,
                        aperiodic_rdegree, periodic_rdegree> const& right);

template<typename LValue, typename RValue,
         int aperiodic_ldegree, int periodic_ldegree,
         int aperiodic_rdegree, int periodic_rdegree>
PiecewisePoissonSeries<Product<LValue, RValue>,
                       std::max({aperiodic_ldegree + aperiodic_rdegree,
                                 aperiodic_ldegree + periodic_rdegree,
                                 periodic_ldegree + aperiodic_rdegree,
                                 periodic_ldegree + periodic_rdegree}),
                       std::max({aperiodic_ldegree + periodic_rdegree,
                                 periodic_ldegree + aperiodic_rdegree,
                                 periodic_ldegree + periodic_rdegree})>
operator*(PoissonSeries<LValue,
                        aperiodic_ldegree, periodic_ldegree> const& left,
          PiecewisePoissonSeries<
              RValue, aperiodic_rdegree, periodic_rdegree> const& right);

template<typename LValue, typename RValue,
         int aperiodic_ldegree, int periodic_ldegree,
         int aperiodic_rdegree, int periodic_rdegree>
PiecewisePoissonSeries<Product<LValue, RValue>,
                       std::max({aperiodic_ldegree + aperiodic_rdegree,
                                 aperiodic_ldegree + periodic_rdegree,
                                 periodic_ldegree + aperiodic_rdegree,
                                 periodic_ldegree + periodic_rdegree}),
                       std::max({aperiodic_ldegree + periodic_rdegree,
                                 periodic_ldegree + aperiodic_rdegree,
                                 periodic_ldegree + periodic_rdegree})>
operator*(PiecewisePoissonSeries<
              LValue, aperiodic_ldegree, periodic_ldegree> const& left,
          PoissonSeries<RValue,
                        aperiodic_rdegree, periodic_rdegree> const& right);

template<
    typename LValue, typename RValue,
    int aperiodic_ldegree, int periodic_ldegree,
    int aperiodic_rdegree, int periodic_rdegree,
    int aperiodic_wdegree, int periodic_wdegree>
typename Hilbert<LValue, RValue>::InnerProductType InnerProduct(
    PoissonSeries<LValue,
                  aperiodic_ldegree, periodic_ldegree> const& left,
    PiecewisePoissonSeries<
        RValue, aperiodic_rdegree, periodic_rdegree> const& right,
    PoissonSeries<double,
                  aperiodic_wdegree, periodic_wdegree> const& weight,
    std::optional<int> max_points = std::nullopt);

template<
    typename LValue, typename RValue,
    int aperiodic_ldegree, int periodic_ldegree,
    int aperiodic_rdegree, int periodic_rdegree,
    int aperiodic_wdegree, int periodic_wdegree>
typename Hilbert<LValue, RValue>::InnerProductType InnerProduct(
    PoissonSeries<LValue,
                  aperiodic_ldegree, periodic_ldegree> const& left,
    PiecewisePoissonSeries<RValue,
                           aperiodic_rdegree, periodic_rdegree> const& right,
    PoissonSeries<double,
                  aperiodic_wdegree, periodic_wdegree> const& weight,
    Instant const& t_min,
    Instant const& t_max,
    std::optional<int> max_points = std::nullopt);

template<
    typename LValue, typename RValue,
    int aperiodic_ldegree, int periodic_ldegree,
    int aperiodic_rdegree, int periodic_rdegree,
    int aperiodic_wdegree, int periodic_wdegree>
typename Hilbert<LValue, RValue>::InnerProductType InnerProduct(
    PiecewisePoissonSeries<LValue,
                           aperiodic_ldegree, periodic_ldegree> const& left,
    PoissonSeries<RValue,
                  aperiodic_rdegree, periodic_rdegree> const& right,
    PoissonSeries<double,
                  aperiodic_wdegree, periodic_wdegree> const& weight,
    std::optional<int> max_points = std::nullopt);

template<
    typename LValue, typename RValue,
    int aperiodic_ldegree, int periodic_ldegree,
    int aperiodic_rdegree, int periodic_rdegree,
    int aperiodic_wdegree, int periodic_wdegree>
typename Hilbert<LValue, RValue>::InnerProductType InnerProduct(
    PiecewisePoissonSeries<LValue,
                           aperiodic_ldegree, periodic_ldegree> const& left,
    PoissonSeries<RValue,
                  aperiodic_rdegree, periodic_rdegree> const& right,
    PoissonSeries<double,
                  aperiodic_wdegree, periodic_wdegree> const& weight,
    Instant const& t_min,
    Instant const& t_max,
    std::optional<int> max_points = std::nullopt);

}  // namespace internal

using internal::PiecewisePoissonSeries;

}  // namespace _piecewise_poisson_series
}  // namespace numerics
}  // namespace principia

#include "numerics/piecewise_poisson_series_body.hpp"
