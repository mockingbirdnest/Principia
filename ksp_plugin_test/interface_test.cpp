
#include "ksp_plugin/interface.hpp"

#include <limits>
#include <string>

#include "astronomy/epoch.hpp"
#include "base/not_null.hpp"
#include "base/pull_serializer.hpp"
#include "base/push_deserializer.hpp"
#include "geometry/named_quantities.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "journal/recorder.hpp"
#include "ksp_plugin/frames.hpp"
#include "ksp_plugin/identification.hpp"
#include "ksp_plugin_test/mock_flight_plan.hpp"
#include "ksp_plugin_test/mock_manœuvre.hpp"
#include "ksp_plugin_test/mock_plugin.hpp"
#include "ksp_plugin_test/mock_vessel.hpp"
#include "physics/mock_dynamic_frame.hpp"
#include "quantities/constants.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/almost_equals.hpp"

namespace principia {
namespace interface {

using astronomy::ModifiedJulianDate;
using base::check_not_null;
using base::make_not_null_unique;
using base::PullSerializer;
using base::PushDeserializer;
using geometry::AngularVelocity;
using geometry::Bivector;
using geometry::Displacement;
using geometry::OrthogonalMap;
using geometry::Rotation;
using geometry::Vector;
using geometry::Velocity;
using ksp_plugin::AliceSun;
using ksp_plugin::Barycentric;
using ksp_plugin::Index;
using ksp_plugin::MakeNavigationManœuvre;
using ksp_plugin::MockFlightPlan;
using ksp_plugin::MockManœuvre;
using ksp_plugin::MockPlugin;
using ksp_plugin::MockVessel;
using ksp_plugin::Navball;
using ksp_plugin::Navigation;
using ksp_plugin::NavigationManœuvre;
using ksp_plugin::Part;
using ksp_plugin::PartId;
using ksp_plugin::World;
using ksp_plugin::WorldSun;
using physics::CoordinateFrameField;
using physics::DegreesOfFreedom;
using physics::DynamicFrame;
using physics::Frenet;
using physics::MassiveBody;
using physics::MockDynamicFrame;
using physics::RelativeDegreesOfFreedom;
using physics::RigidMotion;
using physics::RigidTransformation;
using quantities::GravitationalParameter;
using quantities::Length;
using quantities::Pow;
using quantities::SIUnit;
using quantities::Speed;
using quantities::Time;
using quantities::constants::StandardGravity;
using quantities::si::AstronomicalUnit;
using quantities::si::Day;
using quantities::si::Degree;
using quantities::si::Kilo;
using quantities::si::Metre;
using quantities::si::Newton;
using quantities::si::Second;
using quantities::si::Tonne;
using testing_utilities::AlmostEquals;
using ::testing::AllOf;
using ::testing::ByMove;
using ::testing::DoAll;
using ::testing::ElementsAre;
using ::testing::Eq;
using ::testing::ExitedWithCode;
using ::testing::Invoke;
using ::testing::IsNull;
using ::testing::NotNull;
using ::testing::Pointee;
using ::testing::Property;
using ::testing::Ref;
using ::testing::Return;
using ::testing::ReturnRef;
using ::testing::SetArgPointee;
using ::testing::StrictMock;
using ::testing::_;

char const part_name[] = "Picard's chair";
char const vessel_guid[] = "123-456";
char const vessel_name[] = "NCC-1701-D";

Index const celestial_index = 1;
Index const parent_index = 2;
Index const unused = 666;

PartId const part_id = 42;

double const planetarium_rotation = 10;
double const time = 11;

XYZ parent_position = {4, 5, 6};
XYZ parent_velocity = {7, 8, 9};
QP parent_relative_degrees_of_freedom = {parent_position, parent_velocity};

int const trajectory_size = 10;

ACTION_TEMPLATE(FillUniquePtr,
                // Note the comma between int and k:
                HAS_1_TEMPLATE_PARAMS(int, k),
                AND_1_VALUE_PARAMS(ptr)) {
  std::tr1::get<k>(args)->reset(ptr);
}

MATCHER_P4(BurnMatches, thrust, specific_impulse, initial_time, Δv, "") {
  return arg.thrust == thrust &&
         arg.specific_impulse == specific_impulse &&
         arg.initial_time == initial_time &&
         arg.Δv == Δv;
}

class InterfaceTest : public testing::Test {
 protected:
  static void SetUpTestCase() {
    std::string const test_case_name =
        testing::UnitTest::GetInstance()->current_test_case()->name();
    recorder_ = new journal::Recorder(test_case_name + ".journal.hex");
    journal::Recorder::Activate(recorder_);
  }

  static void TearDownTestCase() {
    journal::Recorder::Deactivate();
  }

  InterfaceTest()
      : plugin_(make_not_null_unique<StrictMock<MockPlugin>>()),
        hexadecimal_simple_plugin_(
            ReadFromHexadecimalFile("simple_plugin.proto.hex")),
        serialized_simple_plugin_(
            ReadFromBinaryFile("simple_plugin.proto.bin")) {}

  static std::string ReadFromBinaryFile(std::string const& filename) {
    std::fstream file =
        std::fstream(SOLUTION_DIR / "ksp_plugin_test" / filename,
                     std::ios::in | std::ios::binary);
    CHECK(file.good());
    std::string binary;
    while (!file.eof()) {
      char c;
      file.get(c);
      binary.append(1, c);
    }
    file.close();
    return binary;
  }

  static std::string ReadFromHexadecimalFile(std::string const& filename) {
    std::fstream file =
        std::fstream(SOLUTION_DIR / "ksp_plugin_test" / filename);
    CHECK(file.good());
    std::string hex;
    while (!file.eof()) {
      std::string line;
      std::getline(file, line);
      for (auto const c : line) {
        if ((c >= '0' && c <= '9') || (c >= 'A' && c <= 'F')) {
          hex.append(1, c);
        }
      }
    }
    file.close();
    return hex;
  }

  not_null<std::unique_ptr<StrictMock<MockPlugin>>> plugin_;
  std::string const hexadecimal_simple_plugin_;
  std::string const serialized_simple_plugin_;
  Instant const t0_;
  static journal::Recorder* recorder_;
};

journal::Recorder* InterfaceTest::recorder_ = nullptr;

using InterfaceDeathTest = InterfaceTest;

// And there is only one thing we say to Death.
TEST_F(InterfaceDeathTest, Errors) {
  Plugin* plugin = nullptr;
  EXPECT_DEATH({
    principia__DeletePlugin(nullptr);
  }, "non NULL");
  EXPECT_DEATH({
    principia__UpdateCelestialHierarchy(plugin, celestial_index, parent_index);
  }, "plugin.*non NULL");
  EXPECT_DEATH({
    principia__UpdateCelestialHierarchy(plugin, celestial_index, parent_index);
  }, "plugin.*non NULL");
  EXPECT_DEATH({
    bool inserted;
    principia__InsertOrKeepVessel(plugin,
                                  vessel_guid,
                                  vessel_name,
                                  parent_index,
                                  /*loaded=*/false,
                                  &inserted);
  }, "plugin.*non NULL");
  EXPECT_DEATH({
    principia__InsertUnloadedPart(plugin,
                                  part_id,
                                  part_name,
                                  vessel_guid,
                                  parent_relative_degrees_of_freedom);
  }, "plugin.*non NULL");
  EXPECT_DEATH({
    principia__VesselFromParent(plugin, celestial_index, vessel_guid);
  }, "plugin.*non NULL");
  EXPECT_DEATH({
    principia__CelestialFromParent(plugin, celestial_index);
  }, "plugin.*non NULL");
  EXPECT_DEATH({
    principia__LogFatal("a fatal error");
  }, "a fatal error");
}

TEST_F(InterfaceTest, InitGoogleLogging1) {
  principia__InitGoogleLogging();
}

TEST_F(InterfaceDeathTest, InitGoogleLogging2) {
  // We use EXPECT_EXIT in this test to avoid interfering with the execution of
  // the other tests.
  int const exit_code = 66;

  EXPECT_EXIT({
    google::ShutdownGoogleLogging();
    principia__InitGoogleLogging();
    exit(exit_code);
  }, ExitedWithCode(exit_code), "");
}

TEST_F(InterfaceDeathTest, ActivateRecorder) {
  EXPECT_DEATH({
    journal::Recorder::Deactivate();
    // Fails because the glog directory doesn't exist.
    principia__ActivateRecorder(true);
  }, "glog.Principia.JOURNAL");
}

TEST_F(InterfaceTest, Log) {
  principia__LogInfo("An info");
  principia__LogWarning("A warning");
  principia__LogError("An error");
}

TEST_F(InterfaceTest, NewPlugin) {
  std::unique_ptr<Plugin> plugin(principia__NewPlugin(
                                     "1 s",
                                     "2 s",
                                     planetarium_rotation));
  EXPECT_THAT(plugin, Not(IsNull()));
}

TEST_F(InterfaceTest, DeletePlugin) {
  Plugin const* plugin = plugin_.release();
  principia__DeletePlugin(&plugin);
  EXPECT_THAT(plugin, IsNull());
}

TEST_F(InterfaceTest, InsertMassiveCelestialAbsoluteCartesian) {
  EXPECT_CALL(
      *plugin_,
      InsertCelestialAbsoluteCartesianConstRef(
          celestial_index,
          std::experimental::make_optional(parent_index),
          DegreesOfFreedom<Barycentric>(
              Barycentric::origin +
              Displacement<Barycentric>(
                  {0 * Metre,
                   23.456e-7 * Kilo(Metre),
                   -1 * AstronomicalUnit}),
              Velocity<Barycentric>(
                  {1 * AstronomicalUnit / Day,
                   1 * Kilo(Metre) / Second,
                   1 * Metre / Second})),
          Pointee(
              AllOf(Property(&MassiveBody::is_oblate, false),
                    Property(&MassiveBody::gravitational_parameter,
                             1.2345e6 * SIUnit<GravitationalParameter>())))));
  BodyParameters const body_parameters = {
      "Brian",
      "1.2345e6  m^3/s^2",
      /*reference_instant=*/0.0,
      /*mean_radius=*/"1 m",
      /*axis_right_ascension=*/"0 deg",
      /*axis_declination=*/"90 deg",
      /*reference_angle=*/"0 deg",
      /*angular_velocity=*/"1 rad/s",
      /*j2=*/nullptr,
      /*reference_radius=*/nullptr};
  principia__InsertCelestialAbsoluteCartesian(plugin_.get(),
                                              celestial_index,
                                              &parent_index,
                                              body_parameters,
                                              "0 m",
                                              "23.456e-7 km",
                                              "-1 au",
                                              "1 au / d",
                                              "  1 km/s",
                                              "1  m / s");
}

TEST_F(InterfaceTest, InsertOblateCelestialAbsoluteCartesian) {
  EXPECT_CALL(
      *plugin_,
      InsertCelestialAbsoluteCartesianConstRef(
          celestial_index,
          std::experimental::make_optional(parent_index),
          DegreesOfFreedom<Barycentric>(
              Barycentric::origin +
              Displacement<Barycentric>(
                  {0 * Metre,
                   23.456e-7 * Kilo(Metre),
                   -1 * AstronomicalUnit}),
              Velocity<Barycentric>(
                  {1 * AstronomicalUnit / Day,
                   1 * Kilo(Metre) / Second,
                   1 * Metre / Second})),
          Pointee(
              AllOf(Property(&MassiveBody::is_oblate, true),
                    Property(&MassiveBody::gravitational_parameter,
                             1.2345e6 *
                                 Pow<3>(Kilo(Metre)) / Pow<2>(Second)),
                    Property(&MassiveBody::mean_radius,
                             666 * Kilo(Metre))))));
  BodyParameters const body_parameters = {"that is called Brian",
                                          "1.2345e6  km^3 / s^2",
                                          999.0,
                                          "666 km",
                                          "42 deg",
                                          u8"8°",
                                          "2 rad",
                                          "0.3 rad / d",
                                          "123e-6",
                                          "1000 km"};
  principia__InsertCelestialAbsoluteCartesian(plugin_.get(),
                                              celestial_index,
                                              &parent_index,
                                              body_parameters,
                                              "0 m",
                                              "23.456e-7 km",
                                              "-1 au",
                                              "1 au / d",
                                              "  1 km/s",
                                              "1  m / s");
}

TEST_F(InterfaceTest, UpdateCelestialHierarchy) {
  EXPECT_CALL(*plugin_,
              UpdateCelestialHierarchy(celestial_index, parent_index));
  principia__UpdateCelestialHierarchy(plugin_.get(),
                                      celestial_index,
                                      parent_index);
}

TEST_F(InterfaceTest, EndInitialization) {
  EXPECT_CALL(*plugin_,
              EndInitialization());
  principia__EndInitialization(plugin_.get());
}

TEST_F(InterfaceTest, InsertOrKeepVessel) {
  bool inserted;
  EXPECT_CALL(*plugin_,
              InsertOrKeepVessel(vessel_guid,
                                 vessel_name,
                                 parent_index,
                                 /*loaded=*/false,
                                 Ref(inserted)));
  EXPECT_CALL(*plugin_, HasVessel(vessel_guid))
      .WillOnce(Return(false))
      .WillOnce(Return(true));
  EXPECT_FALSE(plugin_->HasVessel(vessel_guid));
  principia__InsertOrKeepVessel(plugin_.get(),
                                vessel_guid,
                                vessel_name,
                                parent_index,
                                /*loaded=*/false,
                                &inserted);
  EXPECT_TRUE(plugin_->HasVessel(vessel_guid));
}

TEST_F(InterfaceTest, InsertUnloadedPart) {
  EXPECT_CALL(*plugin_,
              InsertUnloadedPart(
                  part_id,
                  part_name,
                  vessel_guid,
                  RelativeDegreesOfFreedom<AliceSun>(
                      Displacement<AliceSun>(
                          {parent_position.x * SIUnit<Length>(),
                           parent_position.y * SIUnit<Length>(),
                           parent_position.z * SIUnit<Length>()}),
                      Velocity<AliceSun>(
                          {parent_velocity.x * SIUnit<Speed>(),
                           parent_velocity.y * SIUnit<Speed>(),
                           parent_velocity.z * SIUnit<Speed>()}))));
  principia__InsertUnloadedPart(plugin_.get(),
                                part_id,
                                part_name,
                                vessel_guid,
                                parent_relative_degrees_of_freedom);
}

TEST_F(InterfaceTest, AdvanceTime) {
  EXPECT_CALL(*plugin_,
              AdvanceTime(t0_ + time * SIUnit<Time>(),
                          planetarium_rotation * Degree));
  principia__AdvanceTime(plugin_.get(), time, planetarium_rotation);
}

TEST_F(InterfaceTest, ForgetAllHistoriesBefore) {
  EXPECT_CALL(*plugin_,
              ForgetAllHistoriesBefore(t0_ + time * SIUnit<Time>()));
  principia__ForgetAllHistoriesBefore(plugin_.get(), time);
}

TEST_F(InterfaceTest, VesselFromParent) {
  EXPECT_CALL(*plugin_,
              VesselFromParent(celestial_index, vessel_guid))
      .WillOnce(Return(RelativeDegreesOfFreedom<AliceSun>(
                           Displacement<AliceSun>(
                               {parent_position.x * SIUnit<Length>(),
                                parent_position.y * SIUnit<Length>(),
                                parent_position.z * SIUnit<Length>()}),
                           Velocity<AliceSun>(
                               {parent_velocity.x * SIUnit<Speed>(),
                                parent_velocity.y * SIUnit<Speed>(),
                                parent_velocity.z * SIUnit<Speed>()}))));
  QP const result = principia__VesselFromParent(plugin_.get(),
                                                celestial_index,
                                                vessel_guid);
  EXPECT_THAT(result, Eq(parent_relative_degrees_of_freedom));
}

TEST_F(InterfaceTest, CelestialFromParent) {
  EXPECT_CALL(*plugin_,
              CelestialFromParent(celestial_index))
      .WillOnce(Return(RelativeDegreesOfFreedom<AliceSun>(
                           Displacement<AliceSun>(
                               {parent_position.x * SIUnit<Length>(),
                                parent_position.y * SIUnit<Length>(),
                                parent_position.z * SIUnit<Length>()}),
                           Velocity<AliceSun>(
                               {parent_velocity.x * SIUnit<Speed>(),
                                parent_velocity.y * SIUnit<Speed>(),
                                parent_velocity.z * SIUnit<Speed>()}))));
  QP const result = principia__CelestialFromParent(plugin_.get(),
                                                    celestial_index);
  EXPECT_THAT(result, Eq(parent_relative_degrees_of_freedom));
}

TEST_F(InterfaceTest, NewNavigationFrame) {
  StrictMock<MockDynamicFrame<Barycentric, Navigation>>* const
      mock_navigation_frame =
          new StrictMock<MockDynamicFrame<Barycentric, Navigation>>;

  NavigationFrameParameters parameters = {
      serialization::BarycentricRotatingDynamicFrame::kExtensionFieldNumber,
      unused,
      celestial_index,
      parent_index};

  EXPECT_CALL(
      *plugin_,
      FillBarycentricRotatingNavigationFrame(celestial_index, parent_index, _))
      .WillOnce(FillUniquePtr<2>(mock_navigation_frame));
  std::unique_ptr<NavigationFrame> navigation_frame(
      principia__NewNavigationFrame(plugin_.get(), parameters));
  EXPECT_EQ(mock_navigation_frame, navigation_frame.get());

  parameters.extension =
      serialization::BodyCentredNonRotatingDynamicFrame::kExtensionFieldNumber;
  parameters.centre_index = celestial_index;

  EXPECT_CALL(*plugin_,
              FillBodyCentredNonRotatingNavigationFrame(celestial_index, _))
      .WillOnce(FillUniquePtr<1>(mock_navigation_frame));
  navigation_frame.release();
  navigation_frame.reset(
      principia__NewNavigationFrame(plugin_.get(), parameters));
  EXPECT_EQ(mock_navigation_frame, navigation_frame.get());
}

TEST_F(InterfaceTest, SetPlottingFrame) {
  StrictMock<MockDynamicFrame<Barycentric, Navigation>>* const
     mock_navigation_frame =
         new StrictMock<MockDynamicFrame<Barycentric, Navigation>>;
  EXPECT_CALL(*plugin_,
              FillBarycentricRotatingNavigationFrame(celestial_index,
                                                     parent_index,
                                                     _))
      .WillOnce(FillUniquePtr<2>(mock_navigation_frame));
  NavigationFrameParameters parameters = {
      serialization::BarycentricRotatingDynamicFrame::kExtensionFieldNumber,
      unused,
      celestial_index,
      parent_index};
  NavigationFrame* navigation_frame =
      principia__NewNavigationFrame(plugin_.get(), parameters);
  EXPECT_EQ(mock_navigation_frame, navigation_frame);
  EXPECT_CALL(*plugin_, SetPlottingFrameConstRef(Ref(*navigation_frame)));
  principia__SetPlottingFrame(plugin_.get(), &navigation_frame);
  EXPECT_THAT(navigation_frame, IsNull());
  EXPECT_CALL(*plugin_, GetPlottingFrame())
      .WillOnce(Return(mock_navigation_frame));
  EXPECT_EQ(mock_navigation_frame, principia__GetPlottingFrame(plugin_.get()));
}

TEST_F(InterfaceTest, RenderedPrediction) {
  StrictMock<MockDynamicFrame<Barycentric, Navigation>>* const
     mock_navigation_frame =
         new StrictMock<MockDynamicFrame<Barycentric, Navigation>>;
  EXPECT_CALL(*plugin_,
              FillBarycentricRotatingNavigationFrame(celestial_index,
                                                    parent_index,
                                                    _))
      .WillOnce(FillUniquePtr<2>(mock_navigation_frame));
  NavigationFrameParameters parameters = {
      serialization::BarycentricRotatingDynamicFrame::kExtensionFieldNumber,
      unused,
      celestial_index,
      parent_index};
  NavigationFrame* navigation_frame =
      principia__NewNavigationFrame(plugin_.get(), parameters);
  EXPECT_EQ(mock_navigation_frame, navigation_frame);

  EXPECT_CALL(*plugin_, SetPlottingFrameConstRef(Ref(*navigation_frame)));
  principia__SetPlottingFrame(plugin_.get(), &navigation_frame);
  EXPECT_THAT(navigation_frame, IsNull());

  // Construct a test rendered trajectory.
  auto rendered_trajectory = make_not_null_unique<DiscreteTrajectory<World>>();
  Position<World> position =
      World::origin + Displacement<World>({1 * SIUnit<Length>(),
                                           2 * SIUnit<Length>(),
                                           3 * SIUnit<Length>()});
  rendered_trajectory->Append(
      t0_, DegreesOfFreedom<World>(position, Velocity<World>()));
  for (int i = 1; i < trajectory_size; ++i) {
    position += Displacement<World>({10 * SIUnit<Length>(),
                                     20 * SIUnit<Length>(),
                                     30 * SIUnit<Length>()});
  rendered_trajectory->Append(
      t0_ + i * Second, DegreesOfFreedom<World>(position, Velocity<World>()));
  }

  EXPECT_CALL(*plugin_,
              FillRenderedPrediction(
                  vessel_guid,
                  World::origin + Displacement<World>(
                                      {parent_position.x * SIUnit<Length>(),
                                       parent_position.y * SIUnit<Length>(),
                                       parent_position.z * SIUnit<Length>()}),
                  _))
      .WillOnce(FillUniquePtr<2>(rendered_trajectory.release()));
  Iterator* iterator =
      principia__RenderedPrediction(plugin_.get(),
                                    vessel_guid,
                                    parent_position);
  EXPECT_EQ(trajectory_size, principia__IteratorSize(iterator));

  // Traverse it and check that we get the right data.
  for (int i = 0; i < trajectory_size; ++i) {
    EXPECT_FALSE(principia__IteratorAtEnd(iterator));
    XYZ const xyz = principia__IteratorGetXYZ(iterator);
    EXPECT_EQ(1 + 10 * i, xyz.x);
    EXPECT_EQ(2 + 20 * i, xyz.y);
    EXPECT_EQ(3 + 30 * i, xyz.z);
    principia__IteratorIncrement(iterator);
  }
  EXPECT_TRUE(principia__IteratorAtEnd(iterator));

  // Delete it.
  EXPECT_THAT(iterator, Not(IsNull()));
  principia__IteratorDelete(&iterator);
  EXPECT_THAT(iterator, IsNull());
}

TEST_F(InterfaceTest, Iterator) {
  StrictMock<MockDynamicFrame<Barycentric, Navigation>>* const
     mock_navigation_frame =
         new StrictMock<MockDynamicFrame<Barycentric, Navigation>>;
  EXPECT_CALL(*plugin_,
              FillBarycentricRotatingNavigationFrame(celestial_index,
                                                     parent_index,
                                                     _))
      .WillOnce(FillUniquePtr<2>(mock_navigation_frame));
  NavigationFrameParameters parameters = {
      serialization::BarycentricRotatingDynamicFrame::kExtensionFieldNumber,
      unused,
      celestial_index,
      parent_index};
  NavigationFrame* navigation_frame =
      principia__NewNavigationFrame(plugin_.get(), parameters);
  EXPECT_EQ(mock_navigation_frame, navigation_frame);

  EXPECT_CALL(*plugin_, SetPlottingFrameConstRef(Ref(*navigation_frame)));
  principia__SetPlottingFrame(plugin_.get(), &navigation_frame);
  EXPECT_THAT(navigation_frame, IsNull());

  // Construct a test rendered trajectory.
  auto rendered_trajectory = make_not_null_unique<DiscreteTrajectory<World>>();
  Position<World> position =
      World::origin + Displacement<World>({1 * SIUnit<Length>(),
                                           2 * SIUnit<Length>(),
                                           3 * SIUnit<Length>()});
  rendered_trajectory->Append(
      t0_, DegreesOfFreedom<World>(position, Velocity<World>()));
  for (int i = 1; i < trajectory_size; ++i) {
    position += Displacement<World>({10 * SIUnit<Length>(),
                                     20 * SIUnit<Length>(),
                                     30 * SIUnit<Length>()});
  rendered_trajectory->Append(
      t0_ + i * Second, DegreesOfFreedom<World>(position, Velocity<World>()));
  }

  // Construct a LineAndIterator.
  EXPECT_CALL(*plugin_,
              FillRenderedVesselTrajectory(
                  vessel_guid,
                  World::origin + Displacement<World>(
                                      {parent_position.x * SIUnit<Length>(),
                                       parent_position.y * SIUnit<Length>(),
                                       parent_position.z * SIUnit<Length>()}),
                  _))
      .WillOnce(FillUniquePtr<2>(rendered_trajectory.release()));
  Iterator* iterator =
      principia__RenderedVesselTrajectory(plugin_.get(),
                                          vessel_guid,
                                          parent_position);
  EXPECT_EQ(trajectory_size, principia__IteratorSize(iterator));

  // Traverse it and check that we get the right data.
  for (int i = 0; i < trajectory_size; ++i) {
    EXPECT_FALSE(principia__IteratorAtEnd(iterator));
    XYZ const xyz = principia__IteratorGetXYZ(iterator);
    EXPECT_EQ(1 + 10 * i, xyz.x);
    EXPECT_EQ(2 + 20 * i, xyz.y);
    EXPECT_EQ(3 + 30 * i, xyz.z);
    principia__IteratorIncrement(iterator);
  }
  EXPECT_TRUE(principia__IteratorAtEnd(iterator));

  // Delete it.
  EXPECT_THAT(iterator, Not(IsNull()));
  principia__IteratorDelete(&iterator);
  EXPECT_THAT(iterator, IsNull());
}

TEST_F(InterfaceTest, PredictionGettersAndSetters) {
  EXPECT_CALL(*plugin_, SetPredictionLength(42 * Second));
  principia__SetPredictionLength(plugin_.get(), 42);
}

TEST_F(InterfaceTest, NavballOrientation) {
  StrictMock<MockDynamicFrame<Barycentric, Navigation>>* const
     mock_navigation_frame =
         new StrictMock<MockDynamicFrame<Barycentric, Navigation>>;
  EXPECT_CALL(*plugin_,
              FillBarycentricRotatingNavigationFrame(celestial_index,
                                                parent_index,
                                                _))
      .WillOnce(FillUniquePtr<2>(mock_navigation_frame));
  NavigationFrameParameters parameters = {
      serialization::BarycentricRotatingDynamicFrame::kExtensionFieldNumber,
      unused,
      celestial_index,
      parent_index};
  NavigationFrame* navigation_frame =
      principia__NewNavigationFrame(plugin_.get(), parameters);
  EXPECT_EQ(mock_navigation_frame, navigation_frame);

  EXPECT_CALL(*plugin_, SetPlottingFrameConstRef(Ref(*navigation_frame)));
  principia__SetPlottingFrame(plugin_.get(), &navigation_frame);
  EXPECT_THAT(navigation_frame, IsNull());

  Position<World> sun_position =
      World::origin + Displacement<World>(
                          {1 * SIUnit<Length>(),
                           2 * SIUnit<Length>(),
                           3 * SIUnit<Length>()});
  EXPECT_CALL(*plugin_, NavballFrameField(sun_position))
      .WillOnce(Return(
          ByMove(std::make_unique<CoordinateFrameField<World, Navball>>())));
  WXYZ q = principia__NavballOrientation(plugin_.get(),
                                         {1, 2, 3},
                                         {2, 3, 5});
  EXPECT_EQ(q.w, 1);
  EXPECT_EQ(q.x, 0);
  EXPECT_EQ(q.y, 0);
  EXPECT_EQ(q.z, 0);
}

TEST_F(InterfaceTest, Frenet) {
  StrictMock<MockDynamicFrame<Barycentric, Navigation>>* const
     mock_navigation_frame =
         new StrictMock<MockDynamicFrame<Barycentric, Navigation>>;
  EXPECT_CALL(*plugin_,
              FillBarycentricRotatingNavigationFrame(celestial_index,
                                                     parent_index,
                                                     _))
      .WillOnce(FillUniquePtr<2>(mock_navigation_frame));
  NavigationFrameParameters parameters = {
      serialization::BarycentricRotatingDynamicFrame::kExtensionFieldNumber,
      unused,
      celestial_index,
      parent_index};
  NavigationFrame* navigation_frame =
      principia__NewNavigationFrame(plugin_.get(), parameters);
  EXPECT_EQ(mock_navigation_frame, navigation_frame);

  EXPECT_CALL(*plugin_, SetPlottingFrameConstRef(Ref(*navigation_frame)));
  principia__SetPlottingFrame(plugin_.get(), &navigation_frame);
  EXPECT_THAT(navigation_frame, IsNull());

  {
    auto const tangent = Vector<double, World>({4, 5, 6});
    EXPECT_CALL(*plugin_, VesselTangent(vessel_guid)).WillOnce(Return(tangent));
    XYZ t = principia__VesselTangent(plugin_.get(), vessel_guid);
    EXPECT_EQ(t.x, tangent.coordinates().x);
    EXPECT_EQ(t.y, tangent.coordinates().y);
    EXPECT_EQ(t.z, tangent.coordinates().z);
  }
  {
    auto const normal = Vector<double, World>({-13, 7, 5});
    EXPECT_CALL(*plugin_, VesselNormal(vessel_guid)).WillOnce(Return(normal));
    XYZ n = principia__VesselNormal(plugin_.get(), vessel_guid);
    EXPECT_EQ(n.x, normal.coordinates().x);
    EXPECT_EQ(n.y, normal.coordinates().y);
    EXPECT_EQ(n.z, normal.coordinates().z);
  }
  {
    auto const binormal = Vector<double, World>({43, 67, 163});
    EXPECT_CALL(*plugin_, VesselBinormal(vessel_guid))
        .WillOnce(Return(binormal));
    XYZ b = principia__VesselBinormal(plugin_.get(), vessel_guid);
    EXPECT_EQ(b.x, binormal.coordinates().x);
    EXPECT_EQ(b.y, binormal.coordinates().y);
    EXPECT_EQ(b.z, binormal.coordinates().z);
  }
  {
    auto const velocity = Velocity<World>(
        {4 * Metre / Second, 5 * Metre / Second, 6 * Metre / Second});
    EXPECT_CALL(*plugin_, VesselVelocity(vessel_guid))
        .WillOnce(Return(velocity));
    XYZ v = principia__VesselVelocity(plugin_.get(), vessel_guid);
    EXPECT_EQ(v.x, velocity.coordinates().x / (Metre / Second));
    EXPECT_EQ(v.y, velocity.coordinates().y / (Metre / Second));
    EXPECT_EQ(v.z, velocity.coordinates().z / (Metre / Second));
  }
}

TEST_F(InterfaceTest, CurrentTime) {
  Instant const mjd0 = ModifiedJulianDate(0);
  EXPECT_CALL(*plugin_, CurrentTime()).WillOnce(Return(mjd0));
  double const current_time = principia__CurrentTime(plugin_.get());
  EXPECT_THAT(t0_ + current_time * Second, Eq(mjd0));
}

TEST_F(InterfaceTest, Apocalypse) {
  char const* details;
  EXPECT_CALL(*plugin_, HasEncounteredApocalypse(_)).WillOnce(Return(false));
  EXPECT_FALSE(principia__HasEncounteredApocalypse(plugin_.get(), &details));

  constexpr char silly_string[] = "silly";
  EXPECT_CALL(*plugin_, HasEncounteredApocalypse(_))
      .WillOnce(DoAll(SetArgPointee<0>(silly_string), Return(true)));
  EXPECT_TRUE(principia__HasEncounteredApocalypse(plugin_.get(), &details));
  EXPECT_STREQ(silly_string, details);
  principia__DeleteString(&details);
  EXPECT_THAT(details, IsNull());
}

TEST_F(InterfaceTest, SerializePlugin) {
  PullSerializer* serializer = nullptr;
  principia::serialization::Plugin message;
  message.ParseFromString(serialized_simple_plugin_);

  EXPECT_CALL(*plugin_, WriteToMessage(_)).WillOnce(SetArgPointee<0>(message));
  char const* serialization =
      principia__SerializePlugin(plugin_.get(), &serializer);
  EXPECT_STREQ(hexadecimal_simple_plugin_.c_str(), serialization);
  EXPECT_EQ(nullptr, principia__SerializePlugin(plugin_.get(), &serializer));
  principia__DeleteString(&serialization);
  EXPECT_THAT(serialization, IsNull());
}

TEST_F(InterfaceTest, DeserializePlugin) {
  PushDeserializer* deserializer = nullptr;
  Plugin const* plugin = nullptr;
  principia__DeserializePlugin(
          hexadecimal_simple_plugin_.c_str(),
          hexadecimal_simple_plugin_.size(),
          &deserializer,
          &plugin);
  principia__DeserializePlugin(hexadecimal_simple_plugin_.c_str(),
                               0,
                               &deserializer,
                               &plugin);
  EXPECT_THAT(plugin, NotNull());
  principia__DeletePlugin(&plugin);
}

TEST_F(InterfaceDeathTest, SettersAndGetters) {
  // We use EXPECT_EXITs in this test to avoid interfering with the execution of
  // the other tests.
  int const exit_code = 66;
  char const exit_message[] = "Exiting";

  EXPECT_EXIT({
    principia__SetBufferedLogging(100);
    ASSERT_EQ(100, principia__GetBufferedLogging());
    std::cerr << exit_message;
    exit(exit_code);
  }, ExitedWithCode(exit_code), exit_message);

  EXPECT_EXIT({
    principia__SetBufferDuration(101);
    ASSERT_EQ(101, principia__GetBufferDuration());
    std::cerr << exit_message;
    exit(exit_code);
  }, ExitedWithCode(exit_code), exit_message);

  EXPECT_EXIT({
    principia__SetSuppressedLogging(102);
    ASSERT_EQ(102, principia__GetSuppressedLogging());
    std::cerr << exit_message;
    exit(exit_code);
  }, ExitedWithCode(exit_code), exit_message);

  EXPECT_EXIT({
    principia__SetVerboseLogging(103);
    ASSERT_EQ(103, principia__GetVerboseLogging());
    std::cerr << exit_message;
    exit(exit_code);
  }, ExitedWithCode(exit_code), exit_message);

  EXPECT_EXIT({
    principia__SetStderrLogging(2);
    ASSERT_EQ(2, principia__GetStderrLogging());
    std::cerr << exit_message;
    exit(exit_code);
  }, ExitedWithCode(exit_code), exit_message);
}

TEST_F(InterfaceTest, FlightPlan) {
  Burn burn = {/*thrust_in_kilonewtons=*/1,
               /*specific_impulse_in_seconds_g0=*/2,
               /*frame=*/{/*extension=*/6000, /*centre=*/celestial_index},
               /*initial_time=*/3,
               /*delta_v=*/{4, 5, 6}};
  StrictMock<MockVessel> vessel;
  StrictMock<MockFlightPlan> flight_plan;

  EXPECT_CALL(*plugin_, HasVessel(vessel_guid))
      .WillRepeatedly(Return(true));
  EXPECT_CALL(*plugin_, GetVessel(vessel_guid))
      .WillRepeatedly(Return(&vessel));
  EXPECT_CALL(vessel, has_flight_plan())
      .WillRepeatedly(Return(true));
  EXPECT_CALL(vessel, flight_plan())
      .WillRepeatedly(ReturnRef(flight_plan));

  EXPECT_TRUE(principia__FlightPlanExists(plugin_.get(), vessel_guid));

  EXPECT_CALL(*plugin_, CreateFlightPlan(vessel_guid,
                                         Instant() + 30 * Second,
                                         100 * Tonne));
  principia__FlightPlanCreate(plugin_.get(),
                              vessel_guid,
                              /*final_time=*/30,
                              /*mass_in_tonnes=*/100);

  EXPECT_CALL(flight_plan, SetDesiredFinalTime(Instant() + 60 * Second))
      .WillOnce(Return(true));
  EXPECT_TRUE(principia__FlightPlanSetDesiredFinalTime(plugin_.get(),
                                                       vessel_guid,
                                                       60));

  EXPECT_CALL(flight_plan, initial_time())
      .WillOnce(Return(Instant() + 3 * Second));
  EXPECT_EQ(3, principia__FlightPlanGetInitialTime(plugin_.get(), vessel_guid));

  EXPECT_CALL(flight_plan, desired_final_time())
      .WillOnce(Return(Instant() + 4 * Second));
  EXPECT_EQ(4, principia__FlightPlanGetDesiredFinalTime(plugin_.get(),
                                                        vessel_guid));

  EXPECT_CALL(
      flight_plan,
      SetAdaptiveStepParameters(AllOf(
          Property(&Ephemeris<Barycentric>::AdaptiveStepParameters::max_steps,
                   11),
          Property(&Ephemeris<Barycentric>::AdaptiveStepParameters::
                       length_integration_tolerance,
                   22 * Metre),
          Property(&Ephemeris<Barycentric>::AdaptiveStepParameters::
                       speed_integration_tolerance,
                   33 * Metre / Second))))
      .WillOnce(Return(true));
  EXPECT_TRUE(principia__FlightPlanSetAdaptiveStepParameters(
                  plugin_.get(),
                  vessel_guid,
                  {/*max_step=*/11,
                   /*length_integration_tolerance=*/22,
                   /*speed_integration_tolerance=*/33}));

  Ephemeris<Barycentric>::AdaptiveStepParameters adaptive_step_parameters(
      DormandElMikkawyPrince1986RKN434FM<Position<Barycentric>>(),
      /*max_steps=*/111,
      /*length_integration_tolerance=*/222 * Metre,
      /*speed_integration_tolerance=*/333 * Metre / Second);
  EXPECT_CALL(flight_plan, adaptive_step_parameters())
      .WillOnce(ReturnRef(adaptive_step_parameters));
  AdaptiveStepParameters expected_adaptive_step_parameters = {
      /*max_step=*/111,
      /*length_integration_tolerance=*/222,
      /*speed_integration_tolerance=*/333};
  EXPECT_EQ(expected_adaptive_step_parameters,
            principia__FlightPlanGetAdaptiveStepParameters(
                plugin_.get(), vessel_guid));

  EXPECT_CALL(*plugin_,
              FillBodyCentredNonRotatingNavigationFrame(celestial_index, _))
      .WillOnce(FillUniquePtr<1>(
                    new StrictMock<MockDynamicFrame<Barycentric, Navigation>>));
  EXPECT_CALL(flight_plan,
              AppendConstRef(
                  BurnMatches(1 * Kilo(Newton),
                              2 * Second * StandardGravity,
                              Instant() + 3 * Second,
                              Velocity<Frenet<Navigation>>(
                                  {4 * (Metre / Second),
                                   5 * (Metre / Second),
                                   6 * (Metre / Second)}))))
      .WillOnce(Return(true));
  EXPECT_TRUE(principia__FlightPlanAppend(plugin_.get(), vessel_guid, burn));

  EXPECT_CALL(flight_plan, number_of_manœuvres())
      .WillOnce(Return(4));
  EXPECT_EQ(4, principia__FlightPlanNumberOfManoeuvres(plugin_.get(),
                                                       vessel_guid));

  auto const plotting_frame =
      make_not_null_unique<MockDynamicFrame<Barycentric, Navigation>>();
  MockDynamicFrame<Barycentric, Navigation> const* const
      navigation_manœuvre_frame =
          new MockDynamicFrame<Barycentric, Navigation>;
  MockManœuvre<Barycentric, Navigation> navigation_manœuvre(
      10 * Kilo(Newton),
      20 * Tonne,
      30 * Second * StandardGravity,
      Vector<double, Frenet<Navigation>>({1, 1, 1}),
      std::unique_ptr<DynamicFrame<Barycentric, Navigation> const>(
          navigation_manœuvre_frame));
  navigation_manœuvre.set_initial_time(Instant());
  navigation_manœuvre.set_duration(7 * Second);
  auto const barycentric_to_plotting = RigidMotion<Barycentric, Navigation>(
      RigidTransformation<Barycentric, Navigation>(
          Barycentric::origin,
          Navigation::origin,
          OrthogonalMap<Barycentric, Navigation>::Identity()),
      AngularVelocity<Barycentric>(),
      Velocity<Barycentric>());
  EXPECT_CALL(flight_plan, GetManœuvre(3))
      .WillOnce(ReturnRef(navigation_manœuvre));
  EXPECT_CALL(*navigation_manœuvre_frame, WriteToMessage(_))
      .WillOnce(Invoke([](not_null<serialization::DynamicFrame*> message) {
        message->MutableExtension(
            serialization::BodyCentredNonRotatingDynamicFrame::extension);
      }));
  EXPECT_CALL(navigation_manœuvre, InertialDirection())
      .WillOnce(Return(Vector<double, Barycentric>({40, 50, 60})));
  EXPECT_CALL(navigation_manœuvre, FrenetFrame())
      .WillOnce(
          Return(OrthogonalMap<Frenet<Navigation>, Barycentric>::Identity()));
  EXPECT_CALL(*plugin_, CurrentTime()).WillOnce(Return(Instant() - 4 * Second));
  EXPECT_CALL(*plugin_, GetPlottingFrame())
      .WillOnce(Return(plotting_frame.get()));
  EXPECT_CALL(*plotting_frame, ToThisFrameAtTime(Instant()))
      .WillOnce(Return(barycentric_to_plotting));
  EXPECT_CALL(*plotting_frame, FromThisFrameAtTime(Instant() - 4 * Second))
      .WillOnce(Return(barycentric_to_plotting.Inverse()));
  EXPECT_CALL(*plugin_, BarycentricToWorldSun())
      .WillOnce(Return(OrthogonalMap<Barycentric, WorldSun>::Identity()));
  auto const navigation_manoeuvre =
      principia__FlightPlanGetManoeuvre(plugin_.get(),
                                        vessel_guid,
                                        3);
  EXPECT_EQ(10, navigation_manoeuvre.burn.thrust_in_kilonewtons);
  EXPECT_EQ(20, navigation_manoeuvre.initial_mass_in_tonnes);
  EXPECT_THAT(navigation_manoeuvre.burn.specific_impulse_in_seconds_g0,
              AlmostEquals(30, 1));
  EXPECT_EQ(40, navigation_manoeuvre.inertial_direction.x);
  EXPECT_EQ(50, navigation_manoeuvre.inertial_direction.y);
  EXPECT_EQ(60, navigation_manoeuvre.inertial_direction.z);

  EXPECT_CALL(flight_plan, number_of_segments())
      .WillOnce(Return(12));
  EXPECT_EQ(12, principia__FlightPlanNumberOfSegments(plugin_.get(),
                                                      vessel_guid));

  auto rendered_trajectory = make_not_null_unique<DiscreteTrajectory<World>>();
  rendered_trajectory->Append(
      t0_, DegreesOfFreedom<World>(World::origin, Velocity<World>()));
  rendered_trajectory->Append(
      t0_ + 1 * Second,
      DegreesOfFreedom<World>(
          World::origin +
              Displacement<World>({0 * Metre, 1 * Metre, 2 * Metre}),
          Velocity<World>()));
  rendered_trajectory->Append(
      t0_ + 2 * Second,
      DegreesOfFreedom<World>(
          World::origin +
              Displacement<World>({0 * Metre, 2 * Metre, 4 * Metre}),
          Velocity<World>()));
  EXPECT_CALL(flight_plan, GetSegment(3, _, _));
  EXPECT_CALL(*plugin_, FillRenderedTrajectoryFromIterators(_, _, _, _))
      .WillOnce(FillUniquePtr<3>(rendered_trajectory.release()));
  auto* const iterator =
      principia__FlightPlanRenderedSegment(plugin_.get(),
                                           vessel_guid,
                                           {0, 1, 2},
                                           3);
  EXPECT_EQ(XYZ({0, 0, 0}), principia__IteratorGetXYZ(iterator));
  principia__IteratorIncrement(iterator);
  EXPECT_EQ(XYZ({0, 1, 2}), principia__IteratorGetXYZ(iterator));
  principia__IteratorIncrement(iterator);
  EXPECT_EQ(XYZ({0, 2, 4}), principia__IteratorGetXYZ(iterator));

  burn.thrust_in_kilonewtons = 10;
  EXPECT_CALL(*plugin_,
              FillBodyCentredNonRotatingNavigationFrame(celestial_index, _))
      .WillOnce(FillUniquePtr<1>(
                    new StrictMock<MockDynamicFrame<Barycentric, Navigation>>));
  EXPECT_CALL(flight_plan,
              ReplaceLastConstRef(
                  BurnMatches(10 * Kilo(Newton),
                              2 * Second * StandardGravity,
                              Instant() + 3 * Second,
                              Velocity<Frenet<Navigation>>(
                                  {4 * (Metre / Second),
                                   5 * (Metre / Second),
                                   6 * (Metre / Second)}))))
      .WillOnce(Return(true));
  EXPECT_TRUE(principia__FlightPlanReplaceLast(plugin_.get(),
                                               vessel_guid,
                                               burn));

  EXPECT_CALL(flight_plan, RemoveLast());
  principia__FlightPlanRemoveLast(plugin_.get(), vessel_guid);

  EXPECT_CALL(vessel, DeleteFlightPlan());
  principia__FlightPlanDelete(plugin_.get(), vessel_guid);
}

}  // namespace interface
}  // namespace principia
