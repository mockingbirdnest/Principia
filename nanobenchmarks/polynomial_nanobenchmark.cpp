#include "geometry/grassmann.hpp"
#include "geometry/instant.hpp"
#include "geometry/space.hpp"
#include "nanobenchmarks/dependencies.hpp"
#include "nanobenchmarks/nanobenchmark.hpp"  // 🧙 For NANOBENCHMARK_*.
#include "numerics/polynomial_in_monomial_basis.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace nanobenchmarks {

using namespace principia::geometry::_grassmann;
using namespace principia::geometry::_instant;
using namespace principia::geometry::_space;
using namespace principia::nanobenchmarks::_dependencies;
using namespace principia::nanobenchmarks::_nanobenchmark;
using namespace principia::numerics::_polynomial_in_monomial_basis;
using namespace principia::physics::_degrees_of_freedom;
using namespace principia::quantities::_named_quantities;
using namespace principia::quantities::_si;

template<typename Value_>
class PolynomialNanobenchmark : public Nanobenchmark<Value_, Instant> {
 protected:
  using P1A = PolynomialInMonomialBasis<Displacement<World>, Instant, 1>;
  using P2A = PolynomialInMonomialBasis<Displacement<World>, Instant, 2>;
  using P3A = PolynomialInMonomialBasis<Displacement<World>, Instant, 3>;
  using P4A = PolynomialInMonomialBasis<Displacement<World>, Instant, 4>;
  using P5A = PolynomialInMonomialBasis<Displacement<World>, Instant, 5>;
  using P10A = PolynomialInMonomialBasis<Displacement<World>, Instant, 10>;
  using P17A = PolynomialInMonomialBasis<Displacement<World>, Instant, 17>;

  PolynomialNanobenchmark()
      : t0_(Instant() + 0.3 * Second),
        c0_({0 * Metre, 0 * Metre, 1 * Metre}),
        c1_({0 * Metre / Second, 1 * Metre / Second, 0 * Metre / Second}),
        c2_({1 * Metre / Second / Second,
             0 * Metre / Second / Second,
             0 * Metre / Second / Second}),
        c3_({1 * Metre / Second / Second / Second,
             0 * Metre / Second / Second / Second,
             0 * Metre / Second / Second / Second}),
        c4_({1 * Metre / Second / Second / Second / Second,
             0 * Metre / Second / Second / Second / Second,
             0 * Metre / Second / Second / Second / Second}),
        c5_({1 * Metre / Second / Second / Second / Second / Second,
             0 * Metre / Second / Second / Second / Second / Second,
             0 * Metre / Second / Second / Second / Second / Second}),
        p1_({c0_, c1_}, t0_),
        p2_({c0_, c1_, c2_}, t0_),
        p3_({c0_, c1_, c2_, c3_}, t0_),
        p4_({c0_, c1_, c2_, c3_, c4_}, t0_),
        p5_({c0_, c1_, c2_, c3_, c4_, c5_}, t0_),
        p10_(P10A::Coefficients{}, t0_),
        p17_(P17A::Coefficients{}, t0_) {};

  Instant const t0_;
  Displacement<World> const c0_;
  Velocity<World> const c1_;
  Vector<Acceleration, World> const c2_;
  Vector<Jerk, World> const c3_;
  Vector<Snap, World> const c4_;
  Vector<Variation<Snap>, World> const c5_;
  P1A const p1_;
  P2A const p2_;
  P3A const p3_;
  P4A const p4_;
  P5A const p5_;
  P10A const p10_;
  P17A const p17_;
};

using PolynomialDisplacementNanobenchmark =
   PolynomialNanobenchmark<Displacement<World>>;
using PolynomialRelativeDegreesOfFreedomNanobenchmark =
   PolynomialNanobenchmark<RelativeDegreesOfFreedom<World>>;

NANOBENCHMARK_FIXTURE(PolynomialDisplacementNanobenchmark, Value01) {
  return p1_(argument);
}

NANOBENCHMARK_FIXTURE(PolynomialDisplacementNanobenchmark, Value02) {
  return p2_(argument);
}

NANOBENCHMARK_FIXTURE(PolynomialDisplacementNanobenchmark, Value03) {
  return p3_(argument);
}

NANOBENCHMARK_FIXTURE(PolynomialDisplacementNanobenchmark, Value04) {
  return p4_(argument);
}

NANOBENCHMARK_FIXTURE(PolynomialDisplacementNanobenchmark, Value05) {
  return p5_(argument);
}

NANOBENCHMARK_FIXTURE(PolynomialDisplacementNanobenchmark, Value10) {
  return p10_(argument);
}

NANOBENCHMARK_FIXTURE(PolynomialDisplacementNanobenchmark, Value17) {
  return p17_(argument);
}

NANOBENCHMARK_FIXTURE(PolynomialRelativeDegreesOfFreedomNanobenchmark,
                      ValueAndDerivative05) {
  return RelativeDegreesOfFreedom<World>(p5_(argument),
                                         p5_.EvaluateDerivative(argument));
}

NANOBENCHMARK_FIXTURE(PolynomialRelativeDegreesOfFreedomNanobenchmark,
                      ValueAndDerivative10) {
  return RelativeDegreesOfFreedom<World>(p10_(argument),
                                         p10_.EvaluateDerivative(argument));
}

NANOBENCHMARK_FIXTURE(PolynomialRelativeDegreesOfFreedomNanobenchmark,
                      ValueAndDerivative17) {
  return RelativeDegreesOfFreedom<World>(p17_(argument),
                                         p17_.EvaluateDerivative(argument));
}

NANOBENCHMARK_FIXTURE(PolynomialRelativeDegreesOfFreedomNanobenchmark,
                      WithDerivative05) {
  Displacement<World> d;
  Velocity<World> v;
  p5_.EvaluateWithDerivative(argument, d, v);
  return RelativeDegreesOfFreedom<World>(d, v);
}

NANOBENCHMARK_FIXTURE(PolynomialRelativeDegreesOfFreedomNanobenchmark,
                      WithDerivative10) {
  Displacement<World> d;
  Velocity<World> v;
  p10_.EvaluateWithDerivative(argument, d, v);
  return RelativeDegreesOfFreedom<World>(d, v);
}

NANOBENCHMARK_FIXTURE(PolynomialRelativeDegreesOfFreedomNanobenchmark,
                      WithDerivative17) {
  Displacement<World> d;
  Velocity<World> v;
  p17_.EvaluateWithDerivative(argument, d, v);
  return RelativeDegreesOfFreedom<World>(d, v);
}

}  // namespace nanobenchmarks
}  // namespace principia
