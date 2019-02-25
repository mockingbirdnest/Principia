#include "astronomy/standard_product_3.hpp"

#include <fstream>
#include <optional>

#include "absl/strings/numbers.h"
#include "absl/strings/str_cat.h"
#include "astronomy/time_scales.hpp"
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
  auto const read_line = [&file, &filename, &line, &line_number, &location]() {
    std::getline(file, line);
    if (file.fail()) {
      CHECK(file.eof()) << "non-EOF failure after " << location;
      LOG(FATAL) << "unexpected end of file after " << location;
    } else {
      ++line_number;
      location = absl::StrCat(filename.string(), " line ", line_number, ": ", line);
    }
  };

  // The specification uses 1-based column indices, and column ranges with
  // bounds included.
  auto const column = [&line, &location, line_number](int index) -> char {
    CHECK_LT(index - 1, line.size()) << location;
    return line[index - 1];
  };
  auto const columns = [&line, &location, line_number](
                            int first, int last) -> std::string_view {
    CHECK_LT(last - 1, line.size()) << location;
    return std::string_view(&line[first - 1], last - first + 1);
  };

  int epochs;
  int number_of_satellites;

  // Header: #a, #b, #c, or #d record.
  read_line();
  CHECK_EQ(column(1), '#') << location;
  CHECK_GE(column(2), 'a') << location;
  CHECK_LE(column(2), 'd') << location;
  version_ = column(2);
  CHECK(absl::SimpleAtoi(columns(33, 39), &epochs)) << location;

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
      auto const satellite_identifier = columns(column, column + 2);
      auto const full_location =
          absl::StrCat(location, " columns ", column, "-", column + 2);
      if (orbits_.size() != number_of_satellites) {
        CHECK_NE(satellite_identifier, "  0") << full_location;
        CHECK(orbits_.emplace(std::piecewise_construct,
                              std::forward_as_tuple(satellite_identifier),
                              std::forward_as_tuple()).second)
            << "duplicate satellite identifier " << satellite_identifier << ": "
            << full_location;
      } else {
        CHECK_EQ(satellite_identifier, "  0") << full_location;
      }
      read_line();
    }
  }
  if (number_of_satellite_id_records < 5) {
    LOG(FATAL) << u8"at least 5 +␣ records expected: " << location;
  }
  if (version_ < 'd' && number_of_satellite_id_records > 5) {
    LOG(FATAL) << u8"exactly 5 +␣ records expected in SP3-"
               << std::string(version_, 1) << ": " << location;
  }

  // Header: ++ records.
  // Ignore the satellite accuracy exponents.
  for (int i = 0; i < number_of_satellite_id_records; ++i) {
    CHECK_EQ(columns(1, 2), "++") << location;
    read_line();
  }

  // Header: first %c record.
  std::function<Instant(std::string const&)> parse_time;
  CHECK_EQ(columns(1, 2), "%c") << location;
  if (version_ < 'c') {
    parse_time = &ParseGPSTime;
  } else {
    auto const time_system = columns(10, 12);
    if (time_system == "GLO" || time_system == "UTC") {
      parse_time = &ParseUTC;
    } else if (time_system == "TAI") {
      parse_time = &ParseTAI;
    } else if (time_system == "BDT") {
      parse_time = &Parse北斗Time;
    } else if (time_system == "GPS" || time_system == "GAL" ||
               time_system == "TAI" || time_system == "IRN" ||
               time_system == "QZS") {
      parse_time = &ParseGPSTime;
    } else {
      LOG(FATAL) << "unexpected time system identifier " << time_system << ": "
                 << location;
    }
  }

  // Header: second %c record.
  read_line();
  CHECK_EQ(columns(1, 2), "%c") << location;

  // Header: %f records.
  read_line();
  CHECK_EQ(columns(1, 2), "%f") << location;
  read_line();
  CHECK_EQ(columns(1, 2), "%f") << location;

  // Header: %i records.
  read_line();
  CHECK_EQ(columns(1, 2), "%i") << location;
  read_line();
  CHECK_EQ(columns(1, 2), "%i") << location;

  // Header: /* records.
  read_line();
  int number_of_comment_records = 0;
  while (columns(1, 2) == "/*") {
    ++number_of_comment_records;
    read_line();
  }
  if (number_of_comment_records < 4) {
    LOG(FATAL) << "at least 4 /* records expected: " << location;
  }
  if (version_ < 'd' && number_of_comment_records > 5) {
    LOG(FATAL) << "exactly 4 /* records expected in SP3-"
               << std::string(version_, 1) << ": " << location;
  }
}

}  // namespace astronomy
}  // namespace principia
