
#pragma once

#include <map>

#include "geometry/named_quantities.hpp"
#include "numerics/polynomial.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"

namespace principia {
namespace numerics {
namespace internal_poisson_series {

using geometry::Instant;
using quantities::AngularFrequency;
using quantities::Product;
using quantities::Quotient;

template<typename Value, int degree_,
         template<typename, typename, int> class Evaluator>
class PoissonSeries {
 public:
  using Polynomial =
      numerics::PolynomialInMonomialBasis<Value, Instant, degree_, Evaluator>;

  struct AngularFrequencyPolynomials {
    Polynomial sin;
    Polynomial cos;
  };

  PoissonSeries(
      Polynomial const& aperiodic,
      std::map<AngularFrequency, AngularFrequencyPolynomials> const& periodic);

 private:
  Polynomial const aperiodic_;
  std::map<AngularFrequency, AngularFrequencyPolynomials> const periodic_;
};

// Vector space of Poisson series.

template<typename Value, int rdegree_,
         template<typename, typename, int> class Evaluator>
constexpr PoissonSeries<Value, rdegree_, Evaluator>
operator+(PoissonSeries<Value, rdegree_, Evaluator> const& right);

template<typename Value, int rdegree_,
         template<typename, typename, int> class Evaluator>
constexpr PoissonSeries<Value, rdegree_, Evaluator>
operator-(PoissonSeries<Value, rdegree_, Evaluator> const& right);

template<typename Value, int ldegree_, int rdegree_,
         template<typename, typename, int> class Evaluator>
constexpr PoissonSeries<Value, std::max(ldegree_, rdegree_), Evaluator>
operator+(PoissonSeries<Value, ldegree_, Evaluator> const& left,
          PoissonSeries<Value, rdegree_, Evaluator> const& right);

template<typename Value, int ldegree_, int rdegree_,
         template<typename, typename, int> class Evaluator>
constexpr PoissonSeries<Value, std::max(ldegree_, rdegree_), Evaluator>
operator-(PoissonSeries<Value, ldegree_, Evaluator> const& left,
          PoissonSeries<Value, rdegree_, Evaluator> const& right);

template<typename Scalar,
         typename Value, int degree_,
         template<typename, typename, int> class Evaluator>
constexpr PoissonSeries<Product<Scalar, Value>, degree_, Evaluator>
operator*(Scalar const& left,
          PoissonSeries<Value, degree_, Evaluator> const& right);

template<typename Scalar,
         typename Value, int degree_,
         template<typename, typename, int> class Evaluator>
constexpr PoissonSeries<Product<Value, Scalar>, degree_, Evaluator>
operator*(PoissonSeries<Value, degree_, Evaluator> const& left,
          Scalar const& right);

template<typename Scalar,
         typename Value, int degree_,
         template<typename, typename, int> class Evaluator>
constexpr PoissonSeries<Quotient<Value, Scalar>, degree_, Evaluator>
operator/(PoissonSeries<Value, degree_, Evaluator> const& left,
          Scalar const& right);

}  // namespace internal_poisson_series

using internal_poisson_series::PoissonSeries;

}  // namespace numerics
}  // namespace principia

#include "numerics/poisson_series_body.hpp"
