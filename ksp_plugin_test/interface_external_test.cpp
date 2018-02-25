
#include "ksp_plugin/interface.hpp"

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "ksp_plugin_test/mock_planetarium.hpp"
#include "ksp_plugin_test/mock_plugin.hpp"
#include "ksp_plugin_test/mock_renderer.hpp"

namespace principia {
namespace interface {

using ksp_plugin::MockPlanetarium;
using ksp_plugin::MockPlugin;
using ::testing::StrictMock;

class InterfaceExternalTest : public ::testing::Test {
 protected:
  InterfaceExternalTest()
      : plugin_(make_not_null_unique<StrictMock<MockPlugin>>()),
        const_plugin_(plugin_.get()) {}

  not_null<std::unique_ptr<StrictMock<MockPlugin>>> const plugin_;
  StrictMock<MockPlugin> const* const const_plugin_;
  Instant const t0_;
};

}  // namespace interface
}  // namespace principia
