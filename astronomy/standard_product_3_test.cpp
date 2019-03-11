
#include "astronomy/standard_product_3.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

using ::testing::Eq;

namespace principia {
namespace astronomy {

class StandardProduct3Test : public ::testing::Test {};

using StandardProduct3DeathTest = StandardProduct3Test;

// Note: in the following, “conformant” means “conformant enough that our parser
// manages to correctly interpret the file”; we may have missed some
// nonconformance, especially since we ignore many fields.

// Conformant position-only files of versions a through d.
TEST_F(StandardProduct3Test, PositionOnly) {
  // SP3-a file from the European Space Operations Centre (European
  // Space Agency).
  StandardProduct3 sp3a(SOLUTION_DIR / "astronomy" / "sp3_orbits" / "esa11802.eph",
                        StandardProduct3::Dialect::Standard);

  // SP3-b file from the laser operational and analysis centre (Лазерный
  // Операционно-Аналитический Центр) Mission Control Centre (Центр Управления
  // Полётами).
  StandardProduct3 sp3b(SOLUTION_DIR / "astronomy" / "sp3_orbits" / "mcc14000.sp3",
                        StandardProduct3::Dialect::Standard);

  // SP3-c and SP3-d file from the Center for Orbit Determination
  // in Europe (Astronomical Institute of the University of Bern).
  // These files are from the Multi-GNSS EXperiment (MGEX).
  StandardProduct3 sp3c(SOLUTION_DIR / "astronomy" / "sp3_orbits" /
                            "COD0MGXFIN_20181260000_01D_05M_ORB.SP3",
                        StandardProduct3::Dialect::Standard);
  StandardProduct3 sp3d(SOLUTION_DIR / "astronomy" / "sp3_orbits" /
                            "COD0MGXFIN_20183640000_01D_05M_ORB.SP3",
                        StandardProduct3::Dialect::Standard);

  EXPECT_THAT(sp3a.version(), Eq(StandardProduct3::Version::A));
  EXPECT_THAT(sp3b.version(), Eq(StandardProduct3::Version::B));
  EXPECT_THAT(sp3c.version(), Eq(StandardProduct3::Version::C));
  EXPECT_THAT(sp3d.version(), Eq(StandardProduct3::Version::D));
}

// Conformant position and velocity files of versions a and c.
TEST_F(StandardProduct3Test, PositionAndVelocity) {
  // SP3-a file from the National Geospatial-intelligence Agency.
  StandardProduct3 sp3a(
      SOLUTION_DIR / "astronomy" / "sp3_orbits" / "nga20342.eph",
      StandardProduct3::Dialect::Standard);

  // SP3-c file from the Centro di Geodesia Spaziale (Agenzia Spaziale
  // Italiana), an ILRS analysis centre.
  // Orbit for ЭТАЛОН-2.
  StandardProduct3 sp3c(SOLUTION_DIR / "astronomy" / "sp3_orbits" /
                            "asi.orb.etalon2.171209.v70.sp3",
                        StandardProduct3::Dialect::Standard);

  EXPECT_THAT(sp3a.version(), Eq(StandardProduct3::Version::A));
  EXPECT_THAT(sp3c.version(), Eq(StandardProduct3::Version::C));
}

// Test that various nonconformant dialects are nonconformant distinct.

TEST_F(StandardProduct3DeathTest, ILRSANonConformance) {
  // We find no /* records, because they are %/* records in ILRSA.
  EXPECT_DEATH(StandardProduct3(SOLUTION_DIR / "astronomy" / "sp3_orbits" /
                                    "ilrsa.orb.lageos2.160319.v35.sp3",
                                StandardProduct3::Dialect::Standard),
               R"(At least 4 /\* records expected)");
  // We fail to parse the date in the epoch header record, as the fields of
  // ILRSA are correctly aligned, but we expect the ILRSB misalignment.
  // date_time_body.hpp.
  EXPECT_DEATH(
      StandardProduct3(SOLUTION_DIR / "astronomy" / "sp3_orbits" /
                           "ilrsa.orb.lageos2.160319.v35.sp3",
                       StandardProduct3::Dialect::ILRSB),
      R"(SimpleAtoi.* line 23: \*  2016  3 13  0  0  0.00000000 columns 17-18)");
}

TEST_F(StandardProduct3DeathTest, ILRSBNonConformance) {
  // We find no /* records, because they are %/* records in ILRSB.
  EXPECT_DEATH(StandardProduct3(SOLUTION_DIR / "astronomy" / "sp3_orbits" /
                                    "ilrsb.orb.lageos2.160319.v35.sp3",
                                StandardProduct3::Dialect::Standard),
               R"(At least 4 /\* records expected)");
  // We fail to parse the date in the epoch header record, as the fields are
  // misaligned.  The error reporting is not very helpful: we fail deep in
  // date_time_body.hpp.
  EXPECT_DEATH(StandardProduct3(SOLUTION_DIR / "astronomy" / "sp3_orbits" /
                                    "ilrsb.orb.lageos2.160319.v35.sp3",
                                StandardProduct3::Dialect::ILRSA),
               "date_time_body.hpp");
}

}  // namespace astronomy
}  // namespace principia
