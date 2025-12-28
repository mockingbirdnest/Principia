#include "geometry/frame.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/instant.hpp"
#include "geometry/space.hpp"
#include "nanobenchmarks/nanobenchmark.hpp"  // 🧙 For NANOBENCHMARK_*.
#include "numerics/hermite3.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace nanobenchmarks {

using namespace principia::geometry::_frame;
using namespace principia::geometry::_grassmann;
using namespace principia::geometry::_instant;
using namespace principia::geometry::_space;
using namespace principia::nanobenchmarks::_nanobenchmark;
using namespace principia::numerics::_hermite3;
using namespace principia::quantities::_si;

using World = Frame<serialization::Frame::TestTag,
                    Inertial,
                    Handedness::Right,
                    serialization::Frame::TEST>;

class Hermite3Nanobenchmark
    : public Nanobenchmark<Displacement<World>, Instant> {
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

  Instant const t0_;
  Instant const t1_;
  Displacement<World> const d0_;
  Displacement<World> const d1_;
  Velocity<World> const v0_;
  Velocity<World> const v1_;
  H const h_;
};

#if 0
NANOBENCHMARK_FIXTURE(Hermite3Nanobenchmark,
                               Displacement<World>,
                               Instant,
                               ConstructionAndValue) {
  H const h({t0_, t1_ + x * Second}, {d0_, d1_}, {v0_, v1_});
  return ToDouble(h.Evaluate(t0_));
}
#endif

NANOBENCHMARK_FIXTURE(Hermite3Nanobenchmark, Value) {
  return h_.Evaluate(argument);
}

#if 0
NANOBENCHMARK_FIXTURE(Hermite3Nanobenchmark, ValueAndDerivative) {
  Instant const t = t0_ + x * Second;
  return ToDouble(h_.Evaluate(t), h_.EvaluateDerivative(t));
}

NANOBENCHMARK_FIXTURE(Hermite3Nanobenchmark, WithDerivative) {
  Instant const t = t0_ + x * Second;
  Displacement<World> d;
  Velocity<World> v;
  h_.EvaluateWithDerivative(t, d, v);
  return ToDouble(d, v);
}
#endif

using VideNanobenchmark = Nanobenchmark<Displacement<World>, Instant>;

NANOBENCHMARK_EXTERN_ALTERNATE_NAME_FUNCTION(
    VideNanobenchmark,
    googoogoo,
    "?googoogoo@nanobenchmarks@principia@@YQ?AV?$Multivector@V?$Quantity@U?$"
    "Dimensions@$00$0A@$0A@$0A@$0A@$0A@$0A@$0A@@internal@_dimensions@"
    "quantities@principia@@@internal@_quantities@quantities@principia@@U?$"
    "Frame@W4Frame_TestTag@serialization@principia@@$0A@$00$00@2_frame@"
    "geometry@5@$00@internal@_grassmann@geometry@2@V?$Point@V?$Quantity@U?$"
    "Dimensions@$0A@$0A@$00$0A@$0A@$0A@$0A@$0A@@internal@_dimensions@"
    "quantities@principia@@@internal@_quantities@quantities@principia@@@4_"
    "point@62@@Z")

}  // namespace nanobenchmarks
}  // namespace principia
