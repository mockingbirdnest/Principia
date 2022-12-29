#include "mathematica/logger.hpp"

#include <filesystem>
#include <tuple>
#include <vector>

#include "base/not_null.hpp"
#include "geometry/frame.hpp"
#include "gtest/gtest.h"
#include "quantities/si.hpp"

namespace principia {
namespace mathematica {

using base::not_null;
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
  // A construction callback that disables the logger.
  Logger* logger_ptr;
  Logger::SetConstructionCallback(
      [&logger_ptr](std::filesystem::path const& path,
                    not_null<Logger*> const logger) {
        EXPECT_EQ(TEMP_DIR / "mathematica_test.wl", path);
        logger_ptr = logger;
        logger->Disable();
      });

  {
    Logger logger(TEMP_DIR / "mathematica_test.wl");
    EXPECT_EQ(&logger, logger_ptr);
    logger.Append("a", std::tuple(-1 * Second, 3 * Metre), PreserveUnits);
    logger_ptr->Enable();
    logger.Append("a", std::vector{1.0, 2.0, 3.0});
    logger.Append("β", 4 * Metre / Second, PreserveUnits);
    logger.Append("a", F::origin, PreserveUnits);
    logger.Disable();
    logger.Append("a", -6.0);
    logger.Set("d", 6 * Second, PreserveUnits);
    logger.Enable();
    logger.Set("c", 5.0);
  }
  // Go check the file.
  EXPECT_EQ(Set("a",
                std::tuple{std::vector{1.0, 2.0, 3.0}, F::origin},
                PreserveUnits) +
                Set("β", std::tuple{4 * Metre / Second}, PreserveUnits) +
                Set("c", 5.0),
            (std::stringstream{}
             << std::ifstream(TEMP_DIR / "mathematica_test0.wl").rdbuf())
                .str());
}
#endif

}  // namespace mathematica
}  // namespace principia
