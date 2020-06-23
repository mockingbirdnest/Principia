
#pragma once

#include <map>

#include "numerics/polynomial.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"

namespace principia {
namespace numerics {
namespace internal_poisson_series {

using quantities::AngularFrequency;
using quantities::Product;
using quantities::Quotient;
using quantities::Time;

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

 private:
  Polynomial const aperiodic_;
  PolynomialsByAngularFrequency const periodic_;
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

}  // namespace internal_poisson_series

using internal_poisson_series::PoissonSeries;

}  // namespace numerics
}  // namespace principia

#include "numerics/poisson_series_body.hpp"
