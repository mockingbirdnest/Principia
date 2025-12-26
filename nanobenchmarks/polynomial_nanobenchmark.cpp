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

class PolynomialNanobenchmark
    : public Nanobenchmark<Displacement<World>, Instant> {
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

#if 0
  static double ToDouble(Displacement<World> const& displacement) {
    auto const& coordinates = displacement.coordinates();
    return _mm_cvtsd_f64(
        _mm_and_pd(_mm_set_pd(coordinates.x / Metre, coordinates.y / Metre),
                   _mm_set_sd(coordinates.z / Metre)));
  }

  static double ToDouble(Displacement<World> const& displacement,
                         Velocity<World> const& velocity) {
    auto const& d = displacement.coordinates();
    auto const& v = velocity.coordinates();
    return _mm_cvtsd_f64(
        _mm_and_pd(_mm_and_pd(_mm_set_pd(d.x / Metre, d.y / Metre),
                              _mm_set_pd(d.z / Metre, v.x / (Metre / Second))),
                   _mm_set_pd(v.y / (Metre / Second), v.z / (Metre / Second))));
  }
#endif

  Instant PrepareArgument(double x) const override {
    return t0_ + x * Second;
  }

  double ConsumeValue(Displacement<World> value) const override {
    auto const& coordinates = value.coordinates();
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
  P10A const p10_;
  P17A const p17_;
};

NANOBENCHMARK_FIXTURE_TEMPLATE(PolynomialNanobenchmark,
                               Displacement<World>,
                               Instant,
                               Value01) {
  return p1_(argument);
}

NANOBENCHMARK_FIXTURE_TEMPLATE(PolynomialNanobenchmark,
                               Displacement<World>,
                               Instant,
                               Value02) {
  return p2_(argument);
}

NANOBENCHMARK_FIXTURE_TEMPLATE(PolynomialNanobenchmark,
                               Displacement<World>,
                               Instant,
                               Value03) {
  return p3_(argument);
}

NANOBENCHMARK_FIXTURE_TEMPLATE(PolynomialNanobenchmark,
                               Displacement<World>,
                               Instant,
                               Value04) {
  return p4_(argument);
}

NANOBENCHMARK_FIXTURE_TEMPLATE(PolynomialNanobenchmark,
                               Displacement<World>,
                               Instant,
                               Value05) {
  return p5_(argument);
}

NANOBENCHMARK_FIXTURE_TEMPLATE(PolynomialNanobenchmark,
                               Displacement<World>,
                               Instant,
                               Value10) {
  return p10_(argument);
}

NANOBENCHMARK_FIXTURE_TEMPLATE(PolynomialNanobenchmark,
                               Displacement<World>,
                               Instant,
                               Value17) {
  return p17_(argument);
}

#if 0
NANOBENCHMARK_FIXTURE(PolynomialNanobenchmark, ValueAndDerivative05) {
  Instant const t = t0_ + x * Second;
  return ToDouble(p5_(t), p5_.EvaluateDerivative(t));
}

NANOBENCHMARK_FIXTURE(PolynomialNanobenchmark, ValueAndDerivative10) {
  Instant const t = t0_ + x * Second;
  return ToDouble(p10_(t), p10_.EvaluateDerivative(t));
}

NANOBENCHMARK_FIXTURE(PolynomialNanobenchmark, ValueAndDerivative17) {
  Instant const t = t0_ + x * Second;
  return ToDouble(p17_(t), p17_.EvaluateDerivative(t));
}

NANOBENCHMARK_FIXTURE(PolynomialNanobenchmark, WithDerivative05) {
  Instant const t = t0_ + x * Second;
  Displacement<World> d;
  Velocity<World> v;
  p5_.EvaluateWithDerivative(t, d, v);
  return ToDouble(d, v);
}

NANOBENCHMARK_FIXTURE(PolynomialNanobenchmark, WithDerivative10) {
  Instant const t = t0_ + x * Second;
  Displacement<World> d;
  Velocity<World> v;
  p10_.EvaluateWithDerivative(t, d, v);
  return ToDouble(d, v);
}

NANOBENCHMARK_FIXTURE(PolynomialNanobenchmark, WithDerivative17) {
  Instant const t = t0_ + x * Second;
  Displacement<World> d;
  Velocity<World> v;
  p17_.EvaluateWithDerivative(t, d, v);
  return ToDouble(d, v);
}
#endif

}  // namespace nanobenchmarks
}  // namespace principia
