#pragma once

#include <filesystem>
#include <map>
#include <optional>
#include <string>

#include "astronomy/frames.hpp"
#include "geometry/named_quantities.hpp"

namespace principia {
namespace astronomy {
namespace internal_standard_product_3 {

using geometry::Instant;
using geometry::Position;
using geometry::Velocity;

// A representation of data in the extended standard product 3 orbit format.
// Specification:
// - version a: ftp://igs.org/pub/data/format/sp3_docu.txt;
// - version b: ftp://igs.org/pub/data/format/sp3_glon.txt.
// - version c: ftp://igs.org/pub/data/format/sp3c.txt.
// - version d: ftp://igs.org/pub/data/format/sp3d.pdf.
class StandardProduct3 {
 public:
  enum class Dialect {
    Standard,
    ILRSA,  // %/* for comments, missing EOF.
    ILRSB,  // %/* for comments.
  };

  StandardProduct3(std::filesystem::path const& filename, Dialect dialect);

 private:
  struct OrbitPoint {
    Instant time;
    Position<ITRS> position;
    std::optional<Velocity<ITRS>> velocity;
  };

  // REMOVE BEFORE FLIGHT: Satellites should not be strings, especially since
  // SP3-a "  1" is SP3-b or later "G01".
  std::map<std::string, std::vector<OrbitPoint>, std::less<>> orbits_;

  // 'a' through 'd'.
  char version_;

  bool has_velocities_;
};

}  // namespace internal_standard_product_3

using internal_standard_product_3::StandardProduct3;

}  // namespace astronomy
}  // namespace principia
