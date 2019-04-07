
#include "astronomy/standard_product_3.hpp"

#include <algorithm>
#include <limits>
#include <set>
#include <vector>

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "physics/body_surface_dynamic_frame.hpp"
#include "physics/solar_system.hpp"
#include "testing_utilities/componentwise.hpp"
#include "testing_utilities/numerics.hpp"

namespace principia {
namespace astronomy {

using base::dynamic_cast_not_null;
using base::not_null;
using geometry::Position;
using physics::BodySurfaceDynamicFrame;
using physics::ContinuousTrajectory;
using physics::DiscreteTrajectory;
using integrators::EmbeddedExplicitRungeKuttaNyströmIntegrator;
using integrators::SymmetricLinearMultistepIntegrator;
using integrators::methods::DormandالمكاوىPrince1986RKN434FM;
using integrators::methods::QuinlanTremaine1990Order12;
using physics::DegreesOfFreedom;
using physics::Ephemeris;
using physics::RotatingBody;
using physics::SolarSystem;
using quantities::astronomy::JulianYear;
using quantities::si::Deci;
using quantities::si::Metre;
using quantities::si::Milli;
using quantities::si::Second;
using testing_utilities::AbsoluteError;
using testing_utilities::Componentwise;
using ::testing::AllOf;
using ::testing::Each;
using ::testing::ElementsAre;
using ::testing::Eq;
using ::testing::Field;
using ::testing::Lt;
using ::testing::ResultOf;
using ::testing::SizeIs;
using ::testing::UnorderedElementsAre;
using ::testing::ValuesIn;

class StandardProduct3Test : public ::testing::Test {
 protected:
  StandardProduct3Test() {}

  static std::set<StandardProduct3::SatelliteGroup> SatelliteGroups(
      std::vector<StandardProduct3::SatelliteIdentifier> const& identifiers) {
    std::set<StandardProduct3::SatelliteGroup> result;
    for (auto const& id : identifiers) {
      result.insert(id.group);
    }
    return result;
  }

  static constexpr auto group = &StandardProduct3::SatelliteIdentifier::group;
};

using StandardProduct3DeathTest = StandardProduct3Test;

// Note: in the following, “conformant” means “conformant enough that our parser
// manages to correctly interpret the file”; we may have missed some
// nonconformance, especially since we ignore many fields.

// Conformant position-only files of versions a through d.
TEST_F(StandardProduct3Test, PositionOnly) {
  // SP3-a file from the European Space Operations Centre (European
  // Space Agency).
  StandardProduct3 sp3a(
      SOLUTION_DIR / "astronomy" / "standard_product_3" / "esa11802.eph",
      StandardProduct3::Dialect::Standard);

  // SP3-b file from the laser operational and analysis centre (Лазерный
  // Операционно-Аналитический Центр) Mission Control Centre (Центр Управления
  // Полётами).
  StandardProduct3 sp3b(
      SOLUTION_DIR / "astronomy" / "standard_product_3" / "mcc14000.sp3",
      StandardProduct3::Dialect::Standard);

  // SP3-c and SP3-d file from the Center for Orbit Determination
  // in Europe (Astronomical Institute of the University of Bern).
  // These files are from the Multi-GNSS EXperiment (MGEX).
  StandardProduct3 sp3c(SOLUTION_DIR / "astronomy" / "standard_product_3" /
                            "COD0MGXFIN_20181260000_01D_05M_ORB.SP3",
                        StandardProduct3::Dialect::Standard);
  StandardProduct3 sp3d(SOLUTION_DIR / "astronomy" / "standard_product_3" /
                            "COD0MGXFIN_20183640000_01D_05M_ORB.SP3",
                        StandardProduct3::Dialect::Standard);

  EXPECT_THAT(sp3a.version(), Eq(StandardProduct3::Version::A));
  EXPECT_THAT(sp3b.version(), Eq(StandardProduct3::Version::B));
  EXPECT_THAT(sp3c.version(), Eq(StandardProduct3::Version::C));
  EXPECT_THAT(sp3d.version(), Eq(StandardProduct3::Version::D));

  EXPECT_THAT(
      sp3a.satellites(),
      AllOf(SizeIs(26),
            Each(Field(group, Eq(StandardProduct3::SatelliteGroup::GPS)))));
  EXPECT_THAT(
      sp3b.satellites(),
      AllOf(SizeIs(3),
            Each(Field(group, Eq(StandardProduct3::SatelliteGroup::ГЛОНАСС)))));
  EXPECT_THAT(sp3c.satellites(),
              AllOf(SizeIs(81),
                    ResultOf(&SatelliteGroups,
                             UnorderedElementsAre(
                                 StandardProduct3::SatelliteGroup::北斗,
                                 StandardProduct3::SatelliteGroup::Galileo,
                                 StandardProduct3::SatelliteGroup::GPS,
                                 StandardProduct3::SatelliteGroup::準天頂衛星,
                                 StandardProduct3::SatelliteGroup::ГЛОНАСС))));
  EXPECT_THAT(sp3d.satellites(),
              AllOf(SizeIs(91),
                    ResultOf(&SatelliteGroups,
                             UnorderedElementsAre(
                                 StandardProduct3::SatelliteGroup::北斗,
                                 StandardProduct3::SatelliteGroup::Galileo,
                                 StandardProduct3::SatelliteGroup::GPS,
                                 StandardProduct3::SatelliteGroup::準天頂衛星,
                                 StandardProduct3::SatelliteGroup::ГЛОНАСС))));
}

// Conformant position and velocity files of versions a and c.
TEST_F(StandardProduct3Test, PositionAndVelocity) {
  // SP3-a file from the National Geospatial-intelligence Agency.
  StandardProduct3 sp3a(
      SOLUTION_DIR / "astronomy" / "standard_product_3" / "nga20342.eph",
      StandardProduct3::Dialect::Standard);

  // SP3-c file from the Centro di Geodesia Spaziale (Agenzia Spaziale
  // Italiana), an ILRS analysis centre.
  // Orbit for ЭТАЛОН-2.
  StandardProduct3 sp3c(SOLUTION_DIR / "astronomy" / "standard_product_3" /
                            "asi.orb.etalon2.171209.v70.sp3",
                        StandardProduct3::Dialect::Standard);

  EXPECT_THAT(sp3a.version(), Eq(StandardProduct3::Version::A));
  EXPECT_THAT(sp3c.version(), Eq(StandardProduct3::Version::C));

  EXPECT_THAT(
      sp3a.satellites(),
      AllOf(SizeIs(31),
            Each(Field(group, Eq(StandardProduct3::SatelliteGroup::GPS)))));
  EXPECT_THAT(sp3c.satellites(),
              ElementsAre(StandardProduct3::SatelliteIdentifier{
                  StandardProduct3::SatelliteGroup::General, 54}));
}

// Test that the nonstandard dialects are nonconformant and distinct.

TEST_F(StandardProduct3DeathTest, ILRSANonConformance) {
  // We find no /* records, because they are %/* records in ILRSA.
  EXPECT_DEATH(
      StandardProduct3(SOLUTION_DIR / "astronomy" / "standard_product_3" /
                           "ilrsa.orb.lageos2.160319.v35.sp3",
                       StandardProduct3::Dialect::Standard),
      R"(At least 4 /\* records expected)");
  // We fail to parse the date in the epoch header record, as the fields of
  // ILRSA are correctly aligned, but we expect the ILRSB misalignment.
  // date_time_body.hpp.
  EXPECT_DEATH(
      StandardProduct3(SOLUTION_DIR / "astronomy" / "standard_product_3" /
                           "ilrsa.orb.lageos2.160319.v35.sp3",
                       StandardProduct3::Dialect::ILRSB),
      R"(SimpleAtoi.* line 23: \*  2016  3 13  0  0  0.00000000 columns 17-18)");
}

TEST_F(StandardProduct3DeathTest, ILRSBNonConformance) {
  // We find no /* records, because they are %/* records in ILRSB.
  EXPECT_DEATH(
      StandardProduct3(SOLUTION_DIR / "astronomy" / "standard_product_3" /
                           "ilrsb.orb.lageos2.160319.v35.sp3",
                       StandardProduct3::Dialect::Standard),
      R"(At least 4 /\* records expected)");
  // We fail to parse the date in the epoch header record, as the fields are
  // misaligned.  The error reporting is not very helpful: we fail deep in
  // date_time_body.hpp.
  EXPECT_DEATH(
      StandardProduct3(SOLUTION_DIR / "astronomy" / "standard_product_3" /
                           "ilrsb.orb.lageos2.160319.v35.sp3",
                       StandardProduct3::Dialect::ILRSA),
      "date_time_body.hpp");
}

// Test that we successfully parse the nonstandard dialects.

TEST_F(StandardProduct3Test, Dialects) {
  StandardProduct3 ilrsa(SOLUTION_DIR / "astronomy" / "standard_product_3" /
                             "ilrsa.orb.lageos2.160319.v35.sp3",
                         StandardProduct3::Dialect::ILRSA);

  StandardProduct3 ilrsb(SOLUTION_DIR / "astronomy" / "standard_product_3" /
                             "ilrsb.orb.lageos2.160319.v35.sp3",
                         StandardProduct3::Dialect::ILRSB);

  EXPECT_THAT(ilrsa.satellites(),
              ElementsAre(StandardProduct3::SatelliteIdentifier{
                  StandardProduct3::SatelliteGroup::General, 52}));
  EXPECT_THAT(ilrsb.satellites(),
              ElementsAre(StandardProduct3::SatelliteIdentifier{
                  StandardProduct3::SatelliteGroup::General, 52}));
}

#if !defined(_DEBUG)

struct StandardProduct3Args {
  std::filesystem::path filename;
  StandardProduct3::Dialect dialect;
  StandardProduct3::Version version;
  bool has_velocities;
};

std::ostream& operator<<(std::ostream& out, StandardProduct3Args const& args) {
  return out << args.filename << " interpreted as " << args.dialect
             << " (expected to have version " << args.version << ")";
}

class StandardProduct3DynamicsTest
    : public ::testing::TestWithParam<StandardProduct3Args> {
 protected:
  StandardProduct3DynamicsTest()
      : solar_system_([]() {
          SolarSystem<ICRS> solar_system(
              SOLUTION_DIR / "astronomy" / "sol_gravity_model.proto.txt",
              SOLUTION_DIR / "astronomy" /
                  "sol_initial_state_jd_2436116_311504629.proto.txt");
          std::vector<std::string> names = solar_system.names();
          for (auto const& name : names) {
            if (name != "Earth") {
              solar_system.RemoveMassiveBody(name);
            }
          }
          solar_system.LimitOblatenessToDegree("Earth", 2);
          solar_system.LimitOblatenessToZonal("Earth");
          return solar_system;
        }()),
        // We can use a long time step because, in the absence of other bodies,
        // the Earth will go in a straight line, which will be integrated
        // exactly.
        ephemeris_(solar_system_.MakeEphemeris(
            /*accuracy_parameters=*/{/*fitting_tolerance=*/5 * Milli(Metre),
                                     /*geopotential_tolerance=*/0x1p-24},
            Ephemeris<ICRS>::FixedStepParameters(
                SymmetricLinearMultistepIntegrator<QuinlanTremaine1990Order12,
                                                   Position<ICRS>>(),
                /*step=*/1 * JulianYear))),
        earth_(dynamic_cast_not_null<RotatingBody<ICRS> const*>(
            solar_system_.massive_body(*ephemeris_, "Earth"))),
        earth_trajectory_(*ephemeris_->trajectory(earth_)),
        itrs_(ephemeris_.get(), earth_) {}

  // This ephemeris has only one body, an oblate Earth.
  // The lack of third bodies makes it quick to prolong; it is suitable for
  // minimal sanity checking of Earth orbits.
  // Its |t_min| is the time of the launch of Спутник-1, so it can be used for
  // any artificial satellite.
  SolarSystem<ICRS> const solar_system_;
  not_null<std::unique_ptr<Ephemeris<ICRS>>> const ephemeris_;
  not_null<RotatingBody<ICRS> const*> const earth_;
  ContinuousTrajectory<ICRS> const& earth_trajectory_;
  BodySurfaceDynamicFrame<ICRS, ITRS> itrs_;
};

INSTANTIATE_TEST_CASE_P(
    AllVersionsAndDialects,
    StandardProduct3DynamicsTest,
    ValuesIn(std::vector<StandardProduct3Args>{
        {SOLUTION_DIR / "astronomy" / "standard_product_3" / "esa11802.eph",
         StandardProduct3::Dialect::Standard,
         StandardProduct3::Version::A},
        {SOLUTION_DIR / "astronomy" / "standard_product_3" / "mcc14000.sp3",
         StandardProduct3::Dialect::Standard,
         StandardProduct3::Version::B},
        {SOLUTION_DIR / "astronomy" / "standard_product_3" /
             "COD0MGXFIN_20181260000_01D_05M_ORB.SP3",
         StandardProduct3::Dialect::Standard,
         StandardProduct3::Version::C},
        {SOLUTION_DIR / "astronomy" / "standard_product_3" /
             "COD0MGXFIN_20183640000_01D_05M_ORB.SP3",
         StandardProduct3::Dialect::Standard,
         StandardProduct3::Version::D},
        {SOLUTION_DIR / "astronomy" / "standard_product_3" / "nga20342.eph",
         StandardProduct3::Dialect::Standard,
         StandardProduct3::Version::A},
        {SOLUTION_DIR / "astronomy" / "standard_product_3" /
             "ilrsa.orb.lageos2.160319.v35.sp3",
         StandardProduct3::Dialect::ILRSA,
         StandardProduct3::Version::C},
        {SOLUTION_DIR / "astronomy" / "standard_product_3" /
             "ilrsb.orb.lageos2.160319.v35.sp3",
         StandardProduct3::Dialect::ILRSB,
         StandardProduct3::Version::C},
    }));

// This test checks that, for each point of the orbit, its evolution taking into
// account only a simple oblate Earth is close enough to the next point.
TEST_P(StandardProduct3DynamicsTest, PerturbedKeplerian) {
  StandardProduct3 sp3(GetParam().filename, GetParam().dialect);
  EXPECT_THAT(sp3.version(), Eq(GetParam().version));
  for (auto const& satellite : sp3.satellites()) {
    for (not_null<DiscreteTrajectory<ITRS> const*> const arc :
         sp3.orbit(satellite)) {
      auto it = arc->Begin();
      for (int i = 0;; ++i) {
        DiscreteTrajectory<ICRS> integrated_arc;
        ephemeris_->Prolong(it.time());
        integrated_arc.Append(
            it.time(),
            itrs_.FromThisFrameAtTime(it.time())(it.degrees_of_freedom()));
        if (++it == arc->End()) {
          break;
        }
      ephemeris_->FlowWithAdaptiveStep(
            &integrated_arc,
            Ephemeris<ICRS>::NoIntrinsicAcceleration,
            it.time(),
            Ephemeris<ICRS>::AdaptiveStepParameters(
                EmbeddedExplicitRungeKuttaNyströmIntegrator<
                    DormandالمكاوىPrince1986RKN434FM,
                    Position<ICRS>>(),
                std::numeric_limits<std::int64_t>::max(),
                /*length_integration_tolerance=*/1 * Milli(Metre),
                /*speed_integration_tolerance=*/1 * Milli(Metre) / Second),
            /*max_ephemeris_steps=*/std::numeric_limits<std::int64_t>::max(),
            /*last_point_only=*/true);
      DegreesOfFreedom<ICRS> actual =
          integrated_arc.last().degrees_of_freedom();
      DegreesOfFreedom<ICRS> expected =
          itrs_.FromThisFrameAtTime(it.time())(it.degrees_of_freedom());
      EXPECT_THAT(AbsoluteError(expected.position(), actual.position()),
                  Lt(20 * Metre))
          << "orbit of satellite " << satellite << " flowing from point " << i;
      EXPECT_THAT(AbsoluteError(expected.velocity(), actual.velocity()),
                  Lt(1 * Deci(Metre) / Second))
          << "orbit of satellite " << satellite << " flowing from point " << i;
      }
    }
  }
}

#endif

}  // namespace astronomy
}  // namespace principia
