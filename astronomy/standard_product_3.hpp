#pragma once

#include <filesystem>
#include <map>
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
class StandardProduct3 {
 public:
  enum class Version : char {
    A = 'a',  // Specification: ftp://igs.org/pub/data/format/sp3_docu.txt.
    B = 'b',  // Specification: ftp://igs.org/pub/data/format/sp3_glon.txt.
    C = 'c',  // Specification: ftp://igs.org/pub/data/format/sp3c.txt.
    D = 'd',  // Specification: ftp://igs.org/pub/data/format/sp3d.pdf.
  };
  StandardProduct3(std::filesystem::path const& filename);

 private:
  Version version_;

  struct OrbitPoint {
    Instant time;
    Position<ITRS> position;
    Velocity<ITRS> velocity;
  };

  // REMOVE BEFORE FLIGHT: Satellites should not be strings, especially since
  // SP3-a "01" is SP3-b or later "G01".
  std::map<std::string, std::vector<OrbitPoint>> orbits_;
};

}  // namespace internal_standard_product_3

using internal_standard_product_3::StandardProduct3;

}  // namespace astronomy
}  // namespace principia
