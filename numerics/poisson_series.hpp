#pragma once

#include <algorithm>
#include <functional>
#include <map>
#include <ostream>
#include <string>
#include <utility>
#include <vector>

#include "base/macros.hpp"  // 🧙 For forward declarations.
#include "base/not_null.hpp"
#include "geometry/complexification.hpp"
#include "geometry/hilbert.hpp"
#include "geometry/instant.hpp"
#include "geometry/interval.hpp"
#include "numerics/polynomial_in_monomial_basis.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"
#include "serialization/numerics.pb.h"

namespace principia {
namespace numerics {
FORWARD_DECLARE(
    TEMPLATE(typename Value,
             int aperiodic_degree, int periodic_degree) class,
    PoissonSeries,
    FROM(poisson_series));
}  // namespace numerics

namespace mathematica {
FORWARD_DECLARE_FUNCTION(
    TEMPLATE(typename Value,
             int aperiodic_degree, int periodic_degree,
             typename OptionalExpressIn) std::string,
    ToMathematicaBody,
    (numerics::_poisson_series::PoissonSeries<
         Value, aperiodic_degree, periodic_degree> const& series,
     OptionalExpressIn express_in),
     FROM(mathematica));
}  // namespace mathematica

namespace numerics {
namespace _poisson_series {
namespace internal {

using namespace principia::base::_not_null;
using namespace principia::geometry::_complexification;
using namespace principia::geometry::_hilbert;
using namespace principia::geometry::_instant;
using namespace principia::geometry::_interval;
using namespace principia::numerics::_polynomial_in_monomial_basis;
using namespace principia::quantities::_named_quantities;
using namespace principia::quantities::_quantities;

// A Poisson series is the sum of terms of the form:
//   aₙtⁿ      aₙₖ tⁿ sin ωₖ t      aₙₖ tⁿ cos ωₖ t
// Terms of the first kind are called aperiodic, terms of the second and third
// kind are called periodic.  Poisson series form an algebra that is stable by
// derivation and integration.
template<typename Value,
         int aperiodic_degree_, int periodic_degree_>
class PoissonSeries {
 public:
  using AperiodicPolynomial =
      PolynomialInMonomialBasis<Value, Instant, aperiodic_degree_>;
  using PeriodicPolynomial =
      PolynomialInMonomialBasis<Value, Instant, periodic_degree_>;

  // TODO(phl): If we wanted to have Poisson series returning affine values,
  // these polynomials should be changed to return Difference<Value>.
  struct Polynomials {
    PeriodicPolynomial sin;
    PeriodicPolynomial cos;
  };

  using PolynomialsByAngularFrequency =
      std::vector<std::pair<AngularFrequency, Polynomials>>;

  // The `periodic` vector may contain frequencies that are negative or zero, as
  // well as repeated frequencies.  It is not expected to be ordered.
  PoissonSeries(AperiodicPolynomial const& aperiodic,
                PolynomialsByAngularFrequency const& periodic);

  // A Poisson series may be explicitly converted to a higher degree.
  template<int higher_aperiodic_degree, int higher_periodic_degree>
  explicit operator PoissonSeries<Value,
                                  higher_aperiodic_degree,
                                  higher_periodic_degree>() const;

  template<int aperiodic_rdegree, int periodic_rdegree>
  PoissonSeries& operator+=(
      PoissonSeries<Value, aperiodic_rdegree, periodic_rdegree> const& right);
  template<int aperiodic_rdegree, int periodic_rdegree>
  PoissonSeries& operator-=(
      PoissonSeries<Value, aperiodic_rdegree, periodic_rdegree> const& right);

  Instant const& origin() const;
  AngularFrequency max_ω() const;

  Value operator()(Instant const& t) const;

  // Returns a copy of this series adjusted to the given origin.
  PoissonSeries AtOrigin(Instant const& origin) const;

  // The constant term of the result is zero.
  PoissonSeries<Primitive<Value, Time>,
                aperiodic_degree_ + 1, periodic_degree_ + 1>
  Primitive() const;

  quantities::_named_quantities::Primitive<Value, Time> Integrate(
      Instant const& t1,
      Instant const& t2) const;

  template<int aperiodic_wdegree, int periodic_wdegree>
  typename Hilbert<Value>::NormType Norm(
      PoissonSeries<double,
                    aperiodic_wdegree, periodic_wdegree> const& weight,
      Instant const& t_min,
      Instant const& t_max) const;

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

  // Similar to the previous constructor, except that the `periodic` vector is
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

  template<typename V, int ar, int pr>
  PoissonSeries<V, ar, pr>
  friend operator-(PoissonSeries<V, ar, pr> const& right);
  template<typename V, int al, int pl, int ar, int pr>
  PoissonSeries<V, std::max(al, ar), std::max(pl, pr)>
  friend operator+(PoissonSeries<V, al, pl> const& left,
                   PoissonSeries<V, ar, pr> const& right);
  template<typename V, int al, int pl, int ar, int pr>
  PoissonSeries<V, std::max(al, ar), std::max(pl, pr)>
  friend operator-(PoissonSeries<V, al, pl> const& left,
                   PoissonSeries<V, ar, pr> const& right);
  template<typename S, typename V, int ar, int pr>
  PoissonSeries<Product<S, V>, ar, pr>
  friend operator*(S const& left,
                   PoissonSeries<V, ar, pr> const& right);
  template<typename S, typename V, int al, int pl>
  PoissonSeries<Product<V, S>, al, pl>
  friend operator*(PoissonSeries<V, al, pl> const& left,
                   S const& right);
  template<typename S, typename V, int al, int pl>
  PoissonSeries<Quotient<V, S>, al, pl>
  friend operator/(PoissonSeries<V, al, pl> const& left,
                   S const& right);
  template<typename L, typename R, int al, int pl, int ar, int pr, typename P>
  auto friend Multiply(PoissonSeries<L, al, pl> const& left,
                       PoissonSeries<R, ar, pr> const& right,
                       P const& product);
  template<typename V, int ad, int pd>
  friend std::ostream& operator<<(std::ostream& out,
                                  PoissonSeries<V, ad, pd> const& series);
  template<typename L, typename R,
           int al, int pl, int ar, int pr, int aw, int pw>
  friend typename Hilbert<L, R>::InnerProductType InnerProduct(
      PoissonSeries<L, al, pl> const& left,
      PoissonSeries<R, ar, pr> const& right,
      PoissonSeries<double, aw, pw> const& weight,
      Instant const& t_min,
      Instant const& t_max);
  template<typename V, int ad, int pd, typename O>
  friend std::string mathematica::_mathematica::internal::ToMathematicaBody(
      PoissonSeries<V, ad, pd> const& polynomial,
      O express_in);
  template<typename V, int ad, int pd>
  friend class PoissonSeries;
};

// Vector space of Poisson series.

template<typename Value,
         int aperiodic_rdegree, int periodic_rdegree>
PoissonSeries<Value, aperiodic_rdegree, periodic_rdegree>
operator+(PoissonSeries<Value,
                        aperiodic_rdegree, periodic_rdegree> const& right);

template<typename Value,
         int aperiodic_rdegree, int periodic_rdegree>
PoissonSeries<Value, aperiodic_rdegree, periodic_rdegree>
operator-(PoissonSeries<Value,
                        aperiodic_rdegree, periodic_rdegree> const& right);

template<typename Value,
         int aperiodic_ldegree, int periodic_ldegree,
         int aperiodic_rdegree, int periodic_rdegree>
PoissonSeries<Value,
              std::max(aperiodic_ldegree, aperiodic_rdegree),
              std::max(periodic_ldegree, periodic_rdegree)>
operator+(PoissonSeries<Value,
                        aperiodic_ldegree, periodic_ldegree> const& left,
          PoissonSeries<Value,
                        aperiodic_rdegree, periodic_rdegree> const& right);

template<typename Value,
         int aperiodic_ldegree, int periodic_ldegree,
         int aperiodic_rdegree, int periodic_rdegree>
PoissonSeries<Value,
              std::max(aperiodic_ldegree, aperiodic_rdegree),
              std::max(periodic_ldegree, periodic_rdegree)>
operator-(PoissonSeries<Value,
                        aperiodic_ldegree, periodic_ldegree> const& left,
          PoissonSeries<Value,
                        aperiodic_rdegree, periodic_rdegree> const& right);

template<typename Scalar,
         typename Value,
         int aperiodic_rdegree, int periodic_rdegree>
PoissonSeries<Product<Scalar, Value>,
              aperiodic_rdegree, periodic_rdegree>
operator*(Scalar const& left,
          PoissonSeries<Value,
                        aperiodic_rdegree, periodic_rdegree> const& right);

template<typename Scalar,
         typename Value,
         int aperiodic_ldegree, int periodic_ldegree>
PoissonSeries<Product<Value, Scalar>,
              aperiodic_ldegree, periodic_ldegree>
operator*(PoissonSeries<Value,
                        aperiodic_ldegree, periodic_ldegree> const& left,
          Scalar const& right);

template<typename Scalar,
         typename Value,
         int aperiodic_ldegree, int periodic_ldegree>
PoissonSeries<Quotient<Value, Scalar>,
              aperiodic_ldegree, periodic_ldegree>
operator/(PoissonSeries<Value,
                        aperiodic_ldegree, periodic_ldegree> const& left,
          Scalar const& right);

// Algebra of Poisson series.

template<typename LValue, typename RValue,
         int aperiodic_ldegree, int periodic_ldegree,
         int aperiodic_rdegree, int periodic_rdegree>
PoissonSeries<Product<LValue, RValue>,
              std::max({aperiodic_ldegree + aperiodic_rdegree,
                        aperiodic_ldegree + periodic_rdegree,
                        periodic_ldegree + aperiodic_rdegree,
                        periodic_ldegree + periodic_rdegree}),
              std::max({aperiodic_ldegree + periodic_rdegree,
                        periodic_ldegree + aperiodic_rdegree,
                        periodic_ldegree + periodic_rdegree})>
operator*(PoissonSeries<LValue,
                        aperiodic_ldegree, periodic_ldegree> const& left,
          PoissonSeries<RValue,
                        aperiodic_rdegree, periodic_rdegree> const& right);

// Returns a scalar-valued Poisson series obtained by pointwise inner product of
// two vector-valued series.
template<typename LValue, typename RValue,
         int aperiodic_ldegree, int periodic_ldegree,
         int aperiodic_rdegree, int periodic_rdegree>
PoissonSeries<typename Hilbert<LValue, RValue>::InnerProductType,
              std::max({aperiodic_ldegree + aperiodic_rdegree,
                        aperiodic_ldegree + periodic_rdegree,
                        periodic_ldegree + aperiodic_rdegree,
                        periodic_ldegree + periodic_rdegree}),
              std::max({aperiodic_ldegree + periodic_rdegree,
                        periodic_ldegree + aperiodic_rdegree,
                        periodic_ldegree + periodic_rdegree})>
PointwiseInnerProduct(
    PoissonSeries<LValue,
                  aperiodic_ldegree, periodic_ldegree> const& left,
    PoissonSeries<RValue,
                  aperiodic_rdegree, periodic_rdegree> const& right);

// Output.

template<typename Value,
         int aperiodic_degree, int periodic_degree>
std::ostream& operator<<(
    std::ostream& out,
    PoissonSeries<Value,
                  aperiodic_degree, periodic_degree> const& series);

// Inner product space of Poisson series.

// Technically the weight function must be nonnegative for this to be an inner
// product.  Not sure how this works with the flat-top windows, which can be
// negative.  Note that the result is normalized by dividing by (t_max - t_min).
template<typename LValue, typename RValue,
         int aperiodic_ldegree, int periodic_ldegree,
         int aperiodic_rdegree, int periodic_rdegree,
         int aperiodic_wdegree, int periodic_wdegree>
typename Hilbert<LValue, RValue>::InnerProductType
InnerProduct(PoissonSeries<LValue,
                           aperiodic_ldegree, periodic_ldegree> const& left,
             PoissonSeries<RValue,
                           aperiodic_rdegree, periodic_rdegree> const& right,
             PoissonSeries<double,
                           aperiodic_wdegree, periodic_wdegree> const& weight,
             Instant const& t_min,
             Instant const& t_max);

}  // namespace internal

using internal::InnerProduct;
using internal::PoissonSeries;

}  // namespace _poisson_series
}  // namespace numerics
}  // namespace principia

#include "numerics/poisson_series_body.hpp"
