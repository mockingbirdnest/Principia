#include "ksp_plugin/interface.hpp"

#include <array>
#include <memory>

#include "base/not_null.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "ksp_plugin/frames.hpp"
#include "ksp_plugin/plugin.hpp"
#include "ksp_plugin_test/mock_plugin.hpp"
#include "ksp_plugin_test/mock_renderer.hpp"
#include "physics/mock_rigid_reference_frame.hpp"
#include "physics/rigid_reference_frame.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace interface {

using ::testing::ByMove;
using ::testing::Pointer;
using ::testing::Return;
using ::testing::ReturnRef;
using ::testing::StrictMock;
using namespace principia::base::_not_null;
using namespace principia::ksp_plugin::_frames;
using namespace principia::ksp_plugin::_plugin;
using namespace principia::ksp_plugin_test::_mock_plugin;
using namespace principia::ksp_plugin_test::_mock_renderer;
using namespace principia::physics::_mock_rigid_reference_frame;
using namespace principia::physics::_rigid_reference_frame;
using namespace principia::quantities::_si;

namespace {

constexpr char vessel_guid[] = "123-456";

constexpr Index celestial_index = 1;
constexpr Index parent_index = 2;
constexpr Index unused = 666;

}  // namespace

class InterfaceRendererTest : public ::testing::Test {
 protected:
  InterfaceRendererTest()
    : plugin_(make_not_null_unique<StrictMock<MockPlugin>>()),
      const_plugin_(plugin_.get()) {}

  MockRenderer renderer_;
  not_null<std::unique_ptr<StrictMock<MockPlugin>>> const plugin_;
  StrictMock<MockPlugin> const* const const_plugin_;
  Instant const t0_;
};

TEST_F(InterfaceRendererTest, SetPlottingFrame) {
  auto* const mock_navigation_frame =
      new StrictMock<MockRigidReferenceFrame<Barycentric, Navigation>>;
  EXPECT_CALL(*plugin_,
              NewBarycentricRotatingNavigationFrame(celestial_index,
                                                    parent_index))
      .WillOnce(Return(
          ByMove(std::unique_ptr<
                 StrictMock<MockRigidReferenceFrame<Barycentric, Navigation>>>(
              mock_navigation_frame))));
  std::array<int const*, 2> const celestial_index_array = {&celestial_index,
                                                           nullptr};
  std::array<int const*, 2> const parent_index_array = {&parent_index,
                                                        nullptr};
  PlottingFrameParameters const parameters = {
      .extension = serialization::BarycentricRotatingReferenceFrame::
          kExtensionFieldNumber,
      .centre_index = unused,
      .primary_index = celestial_index_array.data(),
      .secondary_index = parent_index_array.data()};
  EXPECT_CALL(*plugin_, renderer()).WillRepeatedly(ReturnRef(renderer_));
  EXPECT_CALL(*const_plugin_, renderer()).WillRepeatedly(ReturnRef(renderer_));
  EXPECT_CALL(renderer_, SetPlottingFrame(Pointer(mock_navigation_frame)));
  principia__SetPlottingFrame(plugin_.get(), parameters);
}

TEST_F(InterfaceRendererTest, Frenet) {
  auto* const mock_navigation_frame =
      new StrictMock<MockRigidReferenceFrame<Barycentric, Navigation>>;
  EXPECT_CALL(*plugin_,
              NewBarycentricRotatingNavigationFrame(celestial_index,
                                                    parent_index))
      .WillOnce(Return(
          ByMove(std::unique_ptr<
                 StrictMock<MockRigidReferenceFrame<Barycentric, Navigation>>>(
              mock_navigation_frame))));
  std::array<int const*, 2> const celestial_index_array = {&celestial_index,
                                                           nullptr};
  std::array<int const*, 2> const parent_index_array = {&parent_index,
                                                        nullptr};
  PlottingFrameParameters const parameters = {
      .extension = serialization::BarycentricRotatingReferenceFrame::
          kExtensionFieldNumber,
      .centre_index = unused,
      .primary_index = celestial_index_array.data(),
      .secondary_index = parent_index_array.data()};

  EXPECT_CALL(*plugin_, renderer()).WillRepeatedly(ReturnRef(renderer_));
  EXPECT_CALL(renderer_, SetPlottingFrame(Pointer(mock_navigation_frame)));
  principia__SetPlottingFrame(plugin_.get(), parameters);

  {
    auto const tangent = Vector<double, World>({4, 5, 6});
    EXPECT_CALL(*plugin_, VesselTangent(vessel_guid)).WillOnce(Return(tangent));
    XYZ const t = principia__VesselTangent(plugin_.get(), vessel_guid);
    EXPECT_EQ(t.x, tangent.coordinates().x);
    EXPECT_EQ(t.y, tangent.coordinates().y);
    EXPECT_EQ(t.z, tangent.coordinates().z);
  }
  {
    auto const normal = Vector<double, World>({-13, 7, 5});
    EXPECT_CALL(*plugin_, VesselNormal(vessel_guid)).WillOnce(Return(normal));
    XYZ const n = principia__VesselNormal(plugin_.get(), vessel_guid);
    EXPECT_EQ(n.x, normal.coordinates().x);
    EXPECT_EQ(n.y, normal.coordinates().y);
    EXPECT_EQ(n.z, normal.coordinates().z);
  }
  {
    auto const binormal = Vector<double, World>({43, 67, 163});
    EXPECT_CALL(*plugin_, VesselBinormal(vessel_guid))
        .WillOnce(Return(binormal));
    XYZ const b = principia__VesselBinormal(plugin_.get(), vessel_guid);
    EXPECT_EQ(b.x, binormal.coordinates().x);
    EXPECT_EQ(b.y, binormal.coordinates().y);
    EXPECT_EQ(b.z, binormal.coordinates().z);
  }
  {
    auto const velocity = Velocity<World>(
        {4 * Metre / Second, 5 * Metre / Second, 6 * Metre / Second});
    EXPECT_CALL(*plugin_, VesselVelocity(vessel_guid))
        .WillOnce(Return(velocity));
    XYZ const v = principia__VesselVelocity(plugin_.get(), vessel_guid);
    EXPECT_EQ(v.x, velocity.coordinates().x / (Metre / Second));
    EXPECT_EQ(v.y, velocity.coordinates().y / (Metre / Second));
    EXPECT_EQ(v.z, velocity.coordinates().z / (Metre / Second));
  }
}

}  // namespace interface
}  // namespace principia
