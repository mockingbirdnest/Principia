#include "numerics/fast_fourier_transform.hpp"

#include "gtest/gtest.h"

namespace principia {
namespace numerics {

TEST(FastFourierTransformTest, Square) {
  using FFT = FastFourierTransform<std::vector<double>, 8>;
  FFT const transform({1, 1, 1, 1, 0, 0, 0, 0});
  for (auto t : transform.transform_) {
  LOG(ERROR)<<t;
  }
}

}  // namespace numerics
}  // namespace principia
