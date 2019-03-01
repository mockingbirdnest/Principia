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
    // SP3 file conformant with the specification.
    Standard,
    // International Laser Ranging Service official combined product ILRSA.
    // Divergences from the specification:
    // - comments start with %/* instead of /*;
    // - missing last line EOF.
    ILRSA,
    // International Laser Ranging Service official backup combined product
    // ILRSB.
    // Divergences from the specification:
    // - the number of epochs given on the first line is off by one;
    // - comments start with %/* instead of /*;
    // - the fields of the epoch record are shifted one column to the left;
    // - the minute field of the epoch record takes the value 60 at the end of
    //   the hour.
    ILRSB,
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
