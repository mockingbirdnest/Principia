#include "astronomy/standard_product_3.hpp"

#include <fstream>
#include <optional>

#include "absl/strings/str_cat.h"
#include "glog/logging.h"

namespace principia {
namespace astronomy {

StandardProduct3::StandardProduct3(
    std::filesystem::path const& filename) {
  std::ifstream file(filename);
  CHECK(file.good()) << filename;
  std::string line;
  std::getline(file, line);
  for (int line_number = 1;
       !file.eof();
       std::getline(file, line), ++line_number) {
    std::string location =
        absl::StrCat(filename.string(), " line ", line_number, ": ", line);
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

    switch (line_number) {
      case 1:
        CHECK_EQ(column(1), '#') << location;
        version_ = static_cast<Version>(column(2));
        switch(version_) {
          case Version::A:
          case Version::B:
          case Version::C:
          case Version::D:
            break;
          default:
            LOG(FATAL) << "Unexpected version " << std::string(version_, 1);
        }
    }
  }
}

}  // namespace astronomy
}  // namespace principia