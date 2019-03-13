#pragma once

#include <filesystem>
#include <map>
#include <optional>
#include <string>
#include <vector>

#include "astronomy/frames.hpp"
#include "geometry/named_quantities.hpp"
#include "physics/discrete_trajectory.hpp"

namespace principia {
namespace astronomy {
namespace internal_standard_product_3 {

using geometry::Instant;
using geometry::Position;
using geometry::Velocity;
using physics::DiscreteTrajectory;

// A representation of data in the extended standard product 3 orbit format.
// Specification:
// - version a: ftp://igs.org/pub/data/format/sp3_docu.txt;
// - version b: ftp://igs.org/pub/data/format/sp3_glon.txt.
// - version c: ftp://igs.org/pub/data/format/sp3c.txt.
// - version d: ftp://igs.org/pub/data/format/sp3d.pdf.
class StandardProduct3 {
 public:
  enum class Version : char {
    A = 'a',
    B = 'b',
    C = 'c',
    D = 'd',
  };

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

  enum class SatelliteGroup : char {
    General = 'L',
    GPS = 'G',
    ГЛОНАСС = 'R',
    Galileo = 'E',
    北斗 = 'C',
    準天頂衛星 = 'J',
    IRNSS = 'I',
  };

  struct SatelliteIdentifier {
    SatelliteGroup group;
    int index;
  };

  struct OrbitPoint {
    Instant time;
    Position<ITRS> position;
    std::optional<Velocity<ITRS>> velocity;
  };

  StandardProduct3(std::filesystem::path const& filename, Dialect dialect);

  // The satellite identifiers in the order in which they appear in the file
  // (that order is the same in the satellite ID records and within each epoch).
  std::vector<SatelliteIdentifier> const& satellites() const;

  DiscreteTrajectory<ITRS> const& orbit(SatelliteIdentifier const& id) const;

  Version version() const;

 private:
  std::vector<SatelliteIdentifier> satellites_;
  std::map<SatelliteIdentifier, DiscreteTrajectory<ITRS>> orbits_;
  Version version_;

  bool has_velocities_;
};

bool operator==(StandardProduct3::SatelliteIdentifier const& left,
                StandardProduct3::SatelliteIdentifier const& right);

bool operator<(StandardProduct3::SatelliteIdentifier const& left,
               StandardProduct3::SatelliteIdentifier const& right);

std::ostream& operator<<(std::ostream& out,
                         StandardProduct3::Version const& version);

std::ostream& operator<<(std::ostream& out,
                         StandardProduct3::SatelliteGroup const& group);

std::ostream& operator<<(std::ostream& out,
                         StandardProduct3::SatelliteIdentifier const& id);

}  // namespace internal_standard_product_3

using internal_standard_product_3::StandardProduct3;

}  // namespace astronomy
}  // namespace principia
