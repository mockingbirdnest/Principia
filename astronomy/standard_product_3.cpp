#include "astronomy/standard_product_3.hpp"

#include <fstream>
#include <optional>

#include "absl/strings/numbers.h"
#include "absl/strings/str_cat.h"
#include "glog/logging.h"

namespace principia {
namespace astronomy {

StandardProduct3::StandardProduct3(
    std::filesystem::path const& filename) {
  std::ifstream file(filename);
  CHECK(file.good()) << filename;
  std::string line;
  std::string location;
  int line_number = 0;
  auto const read_line = [&line, &line_number, &location]() {
    std::getline(file, line);
    if (file.fail()) {
      CHECK(file.eof()) << "non-EOF failure after " << location;
      LOG(FATAL) << "unexpected end of file after " << location;
    } else {
      ++line_number;
      location = absl::StrCat(filename.string(), " line ", line_number, ": ", line);
    }
  };

  // The specification uses 1-based column indices.
  auto const column = [&line, &location, line_number](int index) -> char {
    CHECK_LT(index - 1, line.size()) << location;
    return column[index - 1];
  };
  auto const columns = [&line, &location, line_number](
                            int first, int last) -> std::string_view {
    CHECK_LT(last - 1, line.size()) << location;
    return std::string_view(line.data()[first - 1], last - first + 1);
  };

  int epochs;
  int number_of_satellites;

  // Header: #a, #b, #c, or #d record.
  read_line();
  CHECK_EQ(column(1), '#') << location;
  version_ = static_cast<Version>(column(2));
  switch(version_) {
    case Version::A:
    case Version::B:
    case Version::C:
    case Version::D:
      break;
    default:
      LOG(FATAL) << "unexpected version " << std::string(version_, 1) << "\n" << location;
  }
  CHECK(absl::SimpleAtoi(columns(33, 39), &epochs)) << location;
  break;

  // Header: ## record.
  read_line();
  CHECK_EQ(columns(1, 2), "##") << location;

  // Header: +␣ records.
  read_line();
  CHECK_EQ(columns(1, 2), "+ ") << location;
  CHECK(absl::SimpleAtoi(columns(4, 6), &number_of_satellites)) << location;

  int number_of_satellite_id_records = 0;
  while (columns(1, 2) == "+ ") {
    ++number_of_satellite_id_records;
    for (int column = 10; column <= 58; column += 3) {
      if (orbits_.size() != number_of_satellites) {
        auto const satellite_identifier = columns(column, column + 2);
        CHECK_NE(satellite_identifier, "  0")
            << location << " columns " << column << "-" << column + 2;
        CHECK(orbits_.emplace(satellite_identifier, {}).second)
            << "duplicate satellite identifier " << satellite_identifier << ": "
            << location << " columns " << column << "-" << column + 2;
      } else {
        CHECK_EQ(satellite_identifier, "  0")
            << location << " columns " << column << "-" << column + 2;
      }
      read_line();
    }
  }
  if (number_of_satellite_id_records < 5) {
    LOG(FATAL) << u8"at least 5 +␣ records expected: " << location;
  }
  if (version_ != Version::D && number_of_satellite_id_records > 5) {
    LOG(FATAL) << u8"exactly 5 +␣ records expected in a pre-SP3-d file: "
               << location;
  }

  // Header: ++ records.
  // Ignore the satellite accuracy exponents.
  for (int i = 0; i < number_of_satellite_id_records; ++i) {
    CHECK(columns(1, 2), "++") << location;
    read_line();
  }
  CHECK(columns(1, 2), "") << location;
}

}  // namespace astronomy
}  // namespace principia