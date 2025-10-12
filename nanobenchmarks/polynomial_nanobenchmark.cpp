#include "geometry/frame.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/instant.hpp"
#include "geometry/space.hpp"
#include "nanobenchmarks/nanobenchmark.hpp"  // 🧙 For NANOBENCHMARK_*.
#include "numerics/polynomial_evaluators.hpp"
#include "numerics/polynomial_in_monomial_basis.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace nanobenchmarks {
namespace _examples {

using namespace principia::geometry::_frame;
using namespace principia::geometry::_grassmann;
using namespace principia::geometry::_instant;
using namespace principia::geometry::_space;
using namespace principia::nanobenchmarks::_nanobenchmark;
using namespace principia::numerics::_polynomial_evaluators;
using namespace principia::numerics::_polynomial_in_monomial_basis;
using namespace principia::quantities::_named_quantities;
using namespace principia::quantities::_si;

using World = Frame<serialization::Frame::TestTag,
                    Inertial,
                    Handedness::Right,
                    serialization::Frame::TEST>;

class PolynomialNanobenchmark : public Nanobenchmark {
 protected:
  using P1A = PolynomialInMonomialBasis<Displacement<World>, Instant, 1>;
  using P2A = PolynomialInMonomialBasis<Displacement<World>, Instant, 2>;
  using P3A = PolynomialInMonomialBasis<Displacement<World>, Instant, 3>;
  using P4A = PolynomialInMonomialBasis<Displacement<World>, Instant, 4>;

  using P5A = PolynomialInMonomialBasis<Displacement<World>, Instant, 5>;
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
        p1_({c0_, c1_}, t0_, with_evaluator<Estrin>),
        p2_({c0_, c1_, c2_}, t0_, with_evaluator<Estrin>),
        p3_({c0_, c1_, c2_, c3_}, t0_, with_evaluator<Estrin>),
        p4_({c0_, c1_, c2_, c3_, c4_}, t0_, with_evaluator<Estrin>),
        p5_({c0_, c1_, c2_, c3_, c4_, c5_}, t0_, with_evaluator<Estrin>) {};

  static double ToDouble(Displacement<World> const& displacement) {
    auto const& coordinates = displacement.coordinates();
    return _mm_cvtsd_f64(
        _mm_and_pd(_mm_set_pd(coordinates.x / Metre, coordinates.y / Metre),
                   _mm_set_sd(coordinates.z / Metre)));
  }

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
};

NANOBENCHMARK_FIXTURE(PolynomialNanobenchmark, Degree1) {
  return ToDouble(p1_(t0_ + x * Second));
}

NANOBENCHMARK_FIXTURE(PolynomialNanobenchmark, Degree2) {
  return ToDouble(p2_(t0_ + x * Second));
}

NANOBENCHMARK_FIXTURE(PolynomialNanobenchmark, Degree3) {
  return ToDouble(p3_(t0_ + x * Second));
}

NANOBENCHMARK_FIXTURE(PolynomialNanobenchmark, Degree4) {
  return ToDouble(p4_(t0_ + x * Second));
}

NANOBENCHMARK_FIXTURE(PolynomialNanobenchmark, Degree5) {
  return ToDouble(p5_(t0_ + x * Second));
}

}  // namespace _examples
}  // namespace nanobenchmarks
}  // namespace principia
