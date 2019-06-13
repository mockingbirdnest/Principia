
#include "numerics/elliptic_functions.hpp"

#include <filesystem>
#include <fstream>
#include <utility>

#include "glog/logging.h"
#include "google/protobuf/io/zero_copy_stream_impl.h"
#include "google/protobuf/text_format.h"
#include "gtest/gtest.h"
#include "serialization/numerics.pb.h"
#include "testing_utilities/almost_equals.hpp"
#include "testing_utilities/is_near.hpp"

namespace principia {

using testing_utilities::AlmostEquals;
using testing_utilities::IsNear;

namespace numerics {

class EllipticFunctionsTest : public ::testing::Test {
 protected:
  serialization::TabulatedData ReadTabulatedData(
      std::filesystem::path const& tabulated_data_filename) {
    serialization::TabulatedData tabulated_data;
    std::ifstream tabulated_data_ifstream(tabulated_data_filename);
    CHECK(tabulated_data_ifstream.good());
    google::protobuf::io::IstreamInputStream tabulated_data_zcs(
        &tabulated_data_ifstream);
    CHECK(google::protobuf::TextFormat::Parse(&tabulated_data_zcs,
                                              &tabulated_data));
    CHECK_LT(0, tabulated_data.entry_size());
    return tabulated_data;
  }
};

TEST_F(EllipticFunctionsTest, Xgscd) {
  // These expected values come from the output given by Fukushima at the end of
  // his xgscd.txt.  For our purpose, an unit test with assertions is more
  // useful than eyeballing the output.
  auto xgscd_expected =
      ReadTabulatedData(SOLUTION_DIR / "numerics" / "xgscd.proto.txt");
  constexpr double PI = 3.1415926535897932384626433;
  constexpr double PIHALF = PI * 0.5;
  double dmc, mc, m, du, u, s, c, d;
  int jend, iend;
  jend = 10;
  iend = 8;
  dmc = 1.0 / static_cast<double>(jend);
  std::printf("%10s%10s%25s%25s%25s\n", "m", "u", "s", "c", "d");
  int expected_index = 0;
  for (int j = 1; j <= jend; ++j) {
    mc = static_cast<double>(j) * dmc;
    m = 1.0 - mc;
    du = Elk(mc) / static_cast<double>(iend);
    for (int i = 0; i <= iend * 8; ++i) {
      u = du * static_cast<double>(i);
      Gscd(u, mc, s, c, d);
      std::printf("%10.5f%10.5f%25.15e%25.15e%25.15e\n", m, u, s, c, d);

      auto const& expected_entry = xgscd_expected.entry(expected_index);
      auto const expected_argument_m = expected_entry.argument(0);
      auto const expected_argument_u = expected_entry.argument(1);
      auto const expected_value_s = expected_entry.value(0);
      auto const expected_value_c = expected_entry.value(1);
      auto const expected_value_d = expected_entry.value(2);
      EXPECT_THAT(m, IsNear(expected_argument_m, 1.001));
      EXPECT_THAT(u, IsNear(expected_argument_u, 1.001));
      EXPECT_THAT(s, AlmostEquals(expected_value_s, 0, 2));
      EXPECT_THAT(c, AlmostEquals(expected_value_c, 0, 3));
      EXPECT_THAT(d, AlmostEquals(expected_value_d, 0, 1));
      ++expected_index;
    }
    std::printf("\n");
  }
}

}  // namespace numerics
}  // namespace principia
