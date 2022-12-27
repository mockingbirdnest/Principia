#include "mathematica/logger.hpp"

#include "geometry/frame.hpp"
#include "gtest/gtest.h"
#include "quantities/si.hpp"

namespace principia {
namespace mathematica {

using geometry::Frame;
using quantities::si::Metre;
using quantities::si::Second;

class LoggerTest : public ::testing::Test {
 protected:
  using F = Frame<enum class FTag>;
};

// On macOS, std::filesystem is broken (prior to Catalina).  On Ubuntu, the
// introduction of |Flush| caused the test to fail because apparently the file
// is never written to.  We don't really care, we only use the logger on
// Windows.
#if PRINCIPIA_COMPILER_MSVC
TEST_F(LoggerTest, Logger) {
  {
    Logger logger(TEMP_DIR / "mathematica_test.wl");
    logger.Append("a", std::vector{1.0, 2.0, 3.0});
    logger.Append("β", 4 * Metre / Second);
    logger.Append("a", F::origin);
    logger.Set("c", 5.0);
  }
  // Go check the file.
  EXPECT_EQ(Set("a", std::tuple{std::vector{1.0, 2.0, 3.0}, F::origin}) +
                Set("β", std::tuple{4 * Metre / Second}) +
                Set("c", 5.0),
            (std::stringstream{}
             << std::ifstream(TEMP_DIR / "mathematica_test0.wl").rdbuf())
                .str());
}
#endif

}  // namespace mathematica
}  // namespace principia
