
#include "ksp_plugin/interface.hpp"

#include <limits>
#include <optional>
#include <string>
#include <vector>

#include "astronomy/time_scales.hpp"
#include "base/not_null.hpp"
#include "base/pull_serializer.hpp"
#include "base/push_deserializer.hpp"
#include "base/serialization.hpp"
#include "geometry/named_quantities.hpp"
#include "google/protobuf/text_format.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "journal/recorder.hpp"
#include "ksp_plugin/frames.hpp"
#include "ksp_plugin/identification.hpp"
#include "ksp_plugin_test/mock_manœuvre.hpp"
#include "ksp_plugin_test/mock_plugin.hpp"
#include "ksp_plugin_test/mock_renderer.hpp"
#include "ksp_plugin_test/mock_vessel.hpp"
#include "physics/mock_dynamic_frame.hpp"
#include "quantities/astronomy.hpp"
#include "quantities/constants.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/actions.hpp"
#include "testing_utilities/almost_equals.hpp"
#include "testing_utilities/matchers.hpp"
#include "testing_utilities/serialization.hpp"

namespace principia {
namespace interface {

using astronomy::operator""_TT;
using base::check_not_null;
using base::make_not_null_unique;
using base::ParseFromBytes;
using base::PullSerializer;
using base::PushDeserializer;
using geometry::AngularVelocity;
using geometry::Bivector;
using geometry::Displacement;
using geometry::OrthogonalMap;
using geometry::RigidTransformation;
using geometry::Rotation;
using geometry::Vector;
using geometry::Velocity;
using ksp_plugin::AliceSun;
using ksp_plugin::Barycentric;
using ksp_plugin::Index;
using ksp_plugin::MakeNavigationManœuvre;
using ksp_plugin::MockManœuvre;
using ksp_plugin::MockPlugin;
using ksp_plugin::MockRenderer;
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
using quantities::GravitationalParameter;
using quantities::Length;
using quantities::Pow;
using quantities::SIUnit;
using quantities::Speed;
using quantities::Time;
using quantities::astronomy::AstronomicalUnit;
using quantities::constants::StandardGravity;
using quantities::si::Day;
using quantities::si::Degree;
using quantities::si::Kilo;
using quantities::si::Metre;
using quantities::si::Newton;
using quantities::si::Second;
using quantities::si::Tonne;
using testing_utilities::AlmostEquals;
using testing_utilities::EqualsProto;
using testing_utilities::FillUniquePtr;
using testing_utilities::ReadFromBinaryFile;
using testing_utilities::ReadFromHexadecimalFile;
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
using ::testing::SetArgReferee;
using ::testing::SetArgPointee;
using ::testing::StrictMock;
using ::testing::_;

namespace {

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

}  // namespace

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
        hexadecimal_simple_plugin_(ReadFromHexadecimalFile(
            SOLUTION_DIR / "ksp_plugin_test" / "simple_plugin.proto.hex")),
        serialized_simple_plugin_(ReadFromBinaryFile(
            SOLUTION_DIR / "ksp_plugin_test" / "simple_plugin.proto.bin")) {}

  not_null<std::unique_ptr<StrictMock<MockPlugin>>> plugin_;
  std::string const hexadecimal_simple_plugin_;
  std::vector<std::uint8_t> const serialized_simple_plugin_;
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
  }, "glog.*Principia.*JOURNAL");
}

TEST_F(InterfaceTest, Log) {
  principia__LogInfo("An info");
  principia__LogWarning("A warning");
  principia__LogError("An error");
}

TEST_F(InterfaceTest, NewPlugin) {
  std::unique_ptr<Plugin> plugin(principia__NewPlugin(
                                     "MJD1",
                                     "MJD2",
                                     planetarium_rotation));
  EXPECT_THAT(plugin, Not(IsNull()));
}

TEST_F(InterfaceTest, DeletePlugin) {
  Plugin const* plugin = plugin_.release();
  principia__DeletePlugin(&plugin);
  EXPECT_THAT(plugin, IsNull());
}

TEST_F(InterfaceTest, InsertMassiveCelestialAbsoluteCartesian) {
  serialization::GravityModel::Body gravity_model;
  CHECK(google::protobuf::TextFormat::ParseFromString(
      R"(name                    : "Brian"
         gravitational_parameter : "1.2345e6  m^3/s^2"
         reference_instant       : "JD2451545"
         mean_radius             : "1 m"
         axis_right_ascension    : "0 deg"
         axis_declination        : "90 deg"
         reference_angle         : "0 deg"
         angular_frequency       : "1 rad/s")",
      &gravity_model));
  serialization::InitialState::Cartesian::Body initial_state;
  CHECK(google::protobuf::TextFormat::ParseFromString(
      R"(name : "Brian"
         x    : "0 m"
         y    : "23.456e-7 km"
         z    : "-1 au"
         vx   : "1 au / d"
         vy   : "  1 km/s"
         vz   : "1  m / s")",
      &initial_state));
  EXPECT_CALL(*plugin_,
              InsertCelestialAbsoluteCartesian(
                  celestial_index,
                  std::make_optional(parent_index),
                  EqualsProto(gravity_model),
                  EqualsProto(initial_state)));

  BodyParameters const body_parameters = {
      "Brian",
      "1.2345e6  m^3/s^2",
      /*reference_instant=*/"JD2451545",
      /*mean_radius=*/"1 m",
      /*axis_right_ascension=*/"0 deg",
      /*axis_declination=*/"90 deg",
      /*reference_angle=*/"0 deg",
      /*angular_velocity=*/"1 rad/s",
      /*reference_radius=*/nullptr,
      /*j2=*/nullptr};
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
  serialization::GravityModel::Body gravity_model;
  CHECK(google::protobuf::TextFormat::ParseFromString(
      u8R"(name                    : "that is called Brian"
           gravitational_parameter : "1.2345e6  km^3 / s^2"
           reference_instant       : "JD2452545"
           mean_radius             : "666 km"
           axis_right_ascension    : "42 deg"
           axis_declination        : "8°"
           reference_angle         : "2 rad"
           angular_frequency       : "0.3 rad / d"
           reference_radius        : "1000 km"
           j2                      : 123e-6)",
      &gravity_model));
  serialization::InitialState::Cartesian::Body initial_state;
  CHECK(google::protobuf::TextFormat::ParseFromString(
      R"(name : "that is called Brian"
         x    : "0 m",
         y    : "23.456e-7 km"
         z    : "-1 au"
         vx   : "1 au / d"
         vy   : "  1 km/s"
         vz   : "1  m / s")",
      &initial_state));
  EXPECT_CALL(*plugin_,
              InsertCelestialAbsoluteCartesian(
                  celestial_index,
                  std::make_optional(parent_index),
                  EqualsProto(gravity_model),
                  EqualsProto(initial_state)));

  BodyParameters const body_parameters = {"that is called Brian",
                                          "1.2345e6  km^3 / s^2",
                                          "JD2452545",
                                          "666 km",
                                          "42 deg",
                                          u8"8°",
                                          "2 rad",
                                          "0.3 rad / d",
                                          "1000 km",
                                          "123e-6"};
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

TEST_F(InterfaceTest, InsertGeopotentialCelestialAbsoluteCartesian) {
  serialization::GravityModel::Body gravity_model;
  CHECK(google::protobuf::TextFormat::ParseFromString(
      u8R"(name                    : "that is called Brian"
           gravitational_parameter : "1.2345e6  km^3 / s^2"
           reference_instant       : "JD2452545"
           mean_radius             : "666 km"
           axis_right_ascension    : "42 deg"
           axis_declination        : "8°"
           reference_angle         : "2 rad"
           angular_frequency       : "0.3 rad / d"
           reference_radius        : "1000 km"
           geopotential            : {
             row {
               degree : 2,
               column {order: 0, cos: 123e-6, sin: 456e-6}
             }
             row {
               degree : 3,
               column {order: 0, cos: 123e-7, sin: -456e-7}
             }
           })",
      &gravity_model));
  serialization::InitialState::Cartesian::Body initial_state;
  CHECK(google::protobuf::TextFormat::ParseFromString(
      R"(name : "that is called Brian"
         x    : "0 m",
         y    : "23.456e-7 km"
         z    : "-1 au"
         vx   : "1 au / d"
         vy   : "  1 km/s"
         vz   : "1  m / s")",
      &initial_state));
  EXPECT_CALL(*plugin_,
              InsertCelestialAbsoluteCartesian(
                  celestial_index,
                  std::make_optional(parent_index),
                  EqualsProto(gravity_model),
                  EqualsProto(initial_state)));

  BodyGeopotentialElement j2 = {"2", "0", "123e-6", "456e-6"};
  BodyGeopotentialElement j3 = {"3", "0", "123e-7", "-456e-7"};
  BodyGeopotentialElement j[] = {j2, j3};
  BodyParameters const body_parameters = {"that is called Brian",
                                          "1.2345e6  km^3 / s^2",
                                          "JD2452545",
                                          "666 km",
                                          "42 deg",
                                          u8"8°",
                                          "2 rad",
                                          "0.3 rad / d",
                                          "1000 km",
                                          /*j2=*/nullptr,
                                          j,
                                          2};
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
                                 Ref(inserted)))
      .WillOnce(SetArgReferee<4>(true));
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

  MockRenderer renderer;
  EXPECT_CALL(*plugin_, renderer()).WillRepeatedly(ReturnRef(renderer));

  EXPECT_CALL(
      *plugin_,
      FillBarycentricRotatingNavigationFrame(celestial_index, parent_index, _))
      .WillOnce(FillUniquePtr<2>(mock_navigation_frame));
  principia__SetPlottingFrame(plugin_.get(), parameters);

  parameters.extension =
      serialization::BodyCentredNonRotatingDynamicFrame::kExtensionFieldNumber;
  parameters.centre_index = celestial_index;

  EXPECT_CALL(*plugin_,
              FillBodyCentredNonRotatingNavigationFrame(celestial_index, _))
      .WillOnce(FillUniquePtr<1>(mock_navigation_frame));
  principia__SetPlottingFrame(plugin_.get(), parameters);
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

  MockRenderer renderer;
  EXPECT_CALL(*plugin_, renderer()).WillRepeatedly(ReturnRef(renderer));
  EXPECT_CALL(renderer, SetPlottingFrameConstRef(Ref(*mock_navigation_frame)));
  principia__SetPlottingFrame(plugin_.get(), parameters);

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

TEST_F(InterfaceTest, CurrentTime) {
  Instant const mjd0 = "MJD0"_TT;
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
  auto const message = ParseFromBytes<principia::serialization::Plugin>(
      serialized_simple_plugin_);

  EXPECT_CALL(*plugin_, WriteToMessage(_)).WillOnce(SetArgPointee<0>(message));
  char const* serialization =
      principia__SerializePlugin(plugin_.get(),
                                 &serializer,
                                 /*compressor=*/nullptr,
                                 /*encoder=*/nullptr);
  EXPECT_STREQ(hexadecimal_simple_plugin_.c_str(), serialization);
  EXPECT_EQ(nullptr,
            principia__SerializePlugin(plugin_.get(),
                                       &serializer,
                                       /*compressor=*/nullptr,
                                       /*encoder=*/nullptr));
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
          &plugin,
          /*compressor=*/nullptr,
          /*encoder=*/nullptr);
  principia__DeserializePlugin(hexadecimal_simple_plugin_.c_str(),
                                          0,
                                          &deserializer,
                                          &plugin,
                                          /*compressor=*/nullptr,
                                          /*encoder=*/nullptr);
  EXPECT_THAT(plugin, NotNull());
  principia__DeletePlugin(&plugin);
}

// Use for debugging saves given by users.
TEST_F(InterfaceTest, DISABLED_DeserializePluginDebug) {
  PushDeserializer* deserializer = nullptr;
  Plugin const* plugin = nullptr;
  std::string const hexadecimal_plugin = ReadFromHexadecimalFile(
      R"(P:\Public Mockingbird\Principia\Crashes\1595\persistent.proto.hex)");
  principia__DeserializePlugin(
          hexadecimal_plugin.c_str(),
          hexadecimal_plugin.size(),
          &deserializer,
          &plugin,
          /*compressor=*/nullptr,
          /*encoder=*/nullptr);
  principia__DeserializePlugin(hexadecimal_plugin.c_str(),
                                          0,
                                          &deserializer,
                                          &plugin,
                                          /*compressor=*/nullptr,
                                          /*encoder=*/nullptr);
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

}  // namespace interface
}  // namespace principia
