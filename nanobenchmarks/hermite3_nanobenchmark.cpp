#include "geometry/grassmann.hpp"
#include "geometry/instant.hpp"
#include "geometry/space.hpp"
#include "nanobenchmarks/dependencies.hpp"
#include "nanobenchmarks/nanobenchmark.hpp"  // 🧙 For NANOBENCHMARK_*.
#include "numerics/hermite3.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace nanobenchmarks {

using namespace principia::geometry::_grassmann;
using namespace principia::geometry::_instant;
using namespace principia::geometry::_space;
using namespace principia::nanobenchmarks::_dependencies;
using namespace principia::nanobenchmarks::_nanobenchmark;
using namespace principia::numerics::_hermite3;
using namespace principia::physics::_degrees_of_freedom;
using namespace principia::quantities::_si;

template<typename Value_>
class Hermite3Nanobenchmark : public Nanobenchmark<Value_, Instant> {
 protected:
  using H = Hermite3<Displacement<World>, Instant>;

  Hermite3Nanobenchmark()
      : t0_(Instant() + 0.3 * Second),
        t1_(Instant() + 3 * Second),
        d0_({0 * Metre, 0 * Metre, 1 * Metre}),
        d1_({1 * Metre, 3 * Metre, 1 * Metre}),
        v0_({0 * Metre / Second, 1 * Metre / Second, 0 * Metre / Second}),
        v1_({0 * Metre / Second, 3 * Metre / Second, 3 * Metre / Second}),
        h_({t0_, t1_}, {d0_, d1_}, {v0_, v1_}) {};

  Instant const t0_;
  Instant const t1_;
  Displacement<World> const d0_;
  Displacement<World> const d1_;
  Velocity<World> const v0_;
  Velocity<World> const v1_;
  H const h_;
};

using Hermite3DisplacementNanobenchmark =
    Hermite3Nanobenchmark<Displacement<World>>;
using Hermite3RelativeDegreesOfFreedomNanobenchmark =
    Hermite3Nanobenchmark<RelativeDegreesOfFreedom<World>>;

NANOBENCHMARK_FIXTURE(Hermite3DisplacementNanobenchmark, ConstructionAndValue) {
  H const h({t0_, t1_}, {d0_, d1_}, {v0_, v1_});
  return h.Evaluate(argument);
}

NANOBENCHMARK_FIXTURE(Hermite3DisplacementNanobenchmark, Value) {
  return h_.Evaluate(argument);
}

NANOBENCHMARK_FIXTURE(Hermite3RelativeDegreesOfFreedomNanobenchmark,
                      ValueAndDerivative) {
  return RelativeDegreesOfFreedom<World>(h_.Evaluate(argument),
                                         h_.EvaluateDerivative(argument));
}

NANOBENCHMARK_FIXTURE(Hermite3RelativeDegreesOfFreedomNanobenchmark,
                      WithDerivative) {
  Displacement<World> d;
  Velocity<World> v;
  h_.EvaluateWithDerivative(argument, d, v);
  return RelativeDegreesOfFreedom<World>(d, v);
}

}  // namespace nanobenchmarks
}  // namespace principia
