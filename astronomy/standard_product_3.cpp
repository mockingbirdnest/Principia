#include "astronomy/standard_product_3.hpp"

#include <algorithm>
#include <fstream>
#include <optional>
#include <string>
#include <vector>

#include "absl/strings/numbers.h"
#include "absl/strings/str_cat.h"
#include "astronomy/time_scales.hpp"
#include "base/map_util.hpp"
#include "glog/logging.h"
#include "numerics/finite_difference.hpp"

namespace principia {
namespace astronomy {
namespace internal_standard_product_3 {

using base::FindOrDie;
using base::make_not_null_unique;
using geometry::Displacement;
using numerics::FiniteDifference;
using quantities::NaN;
using quantities::Speed;
using quantities::si::Deci;
using quantities::si::Kilo;
using quantities::si::Metre;
using quantities::si::Second;

// Given a trajectory whose velocities are bad or absent (e.g., NaN), uses
// n-point finite difference formulæ on the positions to produce a trajectory
// with consistent velocities.
template<int n>
not_null<std::unique_ptr<DiscreteTrajectory<ITRS>>> ComputeVelocities(
    DiscreteTrajectory<ITRS> const& arc) {
  auto result = make_not_null_unique<DiscreteTrajectory<ITRS>>();
  CHECK_GE(arc.Size(), n);
  std::array<Instant, n> times;
  std::array<Position<ITRS>, n> positions;
  auto it = arc.begin();
  for (int k = 0; k < n; ++k, ++it) {
    // TODO(egg): we should check when reading the file that the times are
    // equally spaced at the interval declared in columns 25-38 of SP3
    // line two.
    times[k] = it->time;
    positions[k] = it->degrees_of_freedom.position();
  }
  // We use a central difference formula wherever possible, so we keep
  // |offset| at (n - 1) / 2 except at the beginning and end of the arc.
  int offset = 0;
  for (int i = 0; i < arc.Size(); ++i) {
    result->Append(
        times[offset],
        {positions[offset],
         FiniteDifference(
             /*values=*/positions,
             /*step=*/(times[n - 1] - times[0]) / (n - 1),
             offset)});
    // At every iteration, either |offset| advances, or the |positions|
    // window shifts and |it| advances.
    if (offset < (n - 1) / 2 || it == arc.end()) {
      ++offset;
    } else {
      std::move(positions.begin() + 1, positions.end(), positions.begin());
      std::move(times.begin() + 1, times.end(), times.begin());
      times.back() = it->time;
      positions.back()  = it->degrees_of_freedom.position();
      ++it;
    }
  }
  // Note that having the right number of calls to |Append| does not guarantee
  // this, as appending at an existing time merely emits a warning.
  CHECK_EQ(result->Size(), arc.Size());
  return result;
}

StandardProduct3::StandardProduct3(
    std::filesystem::path const& filename,
    StandardProduct3::Dialect const dialect) {
  std::ifstream file(filename);
  CHECK(file.good()) << filename;
  std::optional<std::string> line;

  std::string location;
  int line_number = 0;
  auto const read_line = [&file, &filename, &line, &line_number, &location]() {
    if (!line.has_value()) {
      line.emplace();
    }
    std::getline(file, *line);
    if (file.fail()) {
      CHECK(file.eof()) << "Non-EOF failure after " << location;
      line.reset();
      location = absl::StrCat(filename.string(), " at end of file");
    } else {
      ++line_number;
      location =
          absl::StrCat(filename.string(), " line ", line_number, ": ", *line);
    }
  };

  // The specification uses 1-based column indices, and column ranges with
  // bounds included.
  auto const column = [&line, &location](int const index) {
    CHECK(line.has_value()) << location;
    CHECK_LT(index - 1, line->size()) << location;
    return (*line)[index - 1];
  };
  auto const columns = [&line, &location](int const first, int const last) {
    CHECK(line.has_value()) << location;
    CHECK_LT(last - 1, line->size()) << location;
    CHECK_LE(first, last) << location;
    return std::string_view(&(*line)[first - 1], last - first + 1);
  };
  auto const float_columns = [&columns, &location](int const first,
                                                   int const last) {
    double result;
    CHECK(absl::SimpleAtod(columns(first, last), &result))
        << location << " columns " << first << "-" << last;
    return result;
  };
  auto const integer_columns = [&columns, &location](int const first,
                                                     int const last) {
    int result;
    CHECK(absl::SimpleAtoi(columns(first, last), &result))
        << location << " columns " << first << "-" << last;
    return result;
  };

  int number_of_epochs;
  int number_of_satellites;

  // Header: # record.
  read_line();
  CHECK_EQ(column(1), '#') << location;
  CHECK_GE(Version{column(2)}, Version::A) << location;
  CHECK_LE(Version{column(2)}, Version::D) << location;
  version_ = Version{column(2)};
  CHECK(column(3) == 'P' || column(3) == 'V') << location;
  has_velocities_ = column(3) == 'V';
  number_of_epochs = integer_columns(33, 39);
  if (dialect == Dialect::ILRSB) {
    --number_of_epochs;
  }

  // Header: ## record.
  read_line();
  CHECK_EQ(columns(1, 2), "##") << location;

  // Header: +␣ records.
  read_line();
  CHECK_EQ(columns(1, 2), "+ ") << location;
  number_of_satellites = integer_columns(4, 6);

  int number_of_satellite_id_records = 0;
  while (columns(1, 2) == "+ ") {
    ++number_of_satellite_id_records;
    for (int c = 10; c <= 58; c += 3) {
      auto const full_location =
          absl::StrCat(location, " columns ", c, "-", c + 2);
      if (orbits_.size() != number_of_satellites) {
        SatelliteIdentifier id;
        if (version_ == Version::A) {
          // Satellite IDs are purely numeric (and implicitly GPS) in SP3-a.
          CHECK_EQ(column(c), ' ') << full_location;
          id.group = SatelliteGroup::GPS;
        } else {
          id.group = SatelliteGroup{column(c)};
          switch (id.group) {
            case SatelliteGroup::GPS:
            case SatelliteGroup::ГЛОНАСС:
              // GPS and ГЛОНАСС (G and R) satellite IDs are supported since
              // SP3-b.
              break;
            case SatelliteGroup::General:
            case SatelliteGroup::Galileo:
            case SatelliteGroup::北斗:
            case SatelliteGroup::みちびき:
            case SatelliteGroup::IRNSS:
              CHECK_GE(version_, Version::C) << full_location;
              break;
            default:
              LOG(FATAL) << "Invalid satellite identifier " << id << ": "
                         << full_location;
          }
        }
        id.index = integer_columns(c + 1, c + 2);
        CHECK_GT(id.index, 0) << full_location;
        auto const [it, inserted] =
            orbits_.emplace(std::piecewise_construct,
                            std::forward_as_tuple(id),
                            std::forward_as_tuple());
        CHECK(inserted) << "Duplicate satellite identifier " << id << ": "
                        << full_location;
        it->second.push_back(make_not_null_unique<DiscreteTrajectory<ITRS>>());
        satellites_.push_back(id);
      } else {
        CHECK_EQ(columns(c, c + 2), "  0") << full_location;
      }
    }
    read_line();
  }
  if (number_of_satellite_id_records < 5) {
    LOG(FATAL) << u8"at least 5 +␣ records expected: " << location;
  }
  if (version_ < Version::D && number_of_satellite_id_records > 5) {
    if (dialect == Dialect::ChineseMGEX) {
      CHECK_EQ(number_of_satellite_id_records, 10)
          << u8"exactly 10 +␣ records expected in the " << dialect << ": "
          << location;
    } else {
      CHECK_EQ(number_of_satellite_id_records, 5)
          << u8"exactly 5 +␣ records expected in SP3-" << version_ << ": "
          << location;
    }
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
  if (version_ < Version::C) {
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
      LOG(FATAL) << "Unexpected time system identifier " << time_system << ": "
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
  while ((dialect == Dialect::ILRSA || dialect == Dialect::ILRSB)
             ? columns(1, 3) == "%/*"
             : columns(1, 2) == "/*") {
    ++number_of_comment_records;
    read_line();
  }
  if (number_of_comment_records < 4) {
    LOG(FATAL) << "At least 4 /* records expected: " << location;
  }
  if (version_ < Version::D && number_of_comment_records > 5) {
    LOG(FATAL) << "Exactly 4 /* records expected in SP3-"
               << version_ << ": " << location;
  }

  for (int i = 0; i < number_of_epochs; ++i) {
    // *␣ record: the epoch header record.
    CHECK_EQ(columns(1, 2), "* ") << location;
    std::string epoch_string;
    if (dialect == Dialect::ILRSB) {
      int minutes = integer_columns(17, 18);
      int hours = integer_columns(14, 15);
      if (minutes == 60) {
        minutes = 0;
        ++hours;
      }
      epoch_string = absl::StrCat(
          columns(3, 6), "-", columns(8, 9), "-", columns(11, 12),
          "T", absl::Dec(hours, absl::kZeroPad2), ":",
          absl::Dec(minutes, absl::kZeroPad2), ":", columns(20, 25));
    } else {
      // Note: the seconds field is an F11.8, spanning columns 21..31, but our
      // time parser only supports milliseconds.
      epoch_string = absl::StrCat(
          columns(4, 7), "-", columns(9, 10), "-", columns(12, 13),
          "T", columns(15, 16), ":", columns(18, 19), ":", columns(21, 26));
    }
    for (char& c : epoch_string) {
      if (c == ' ') {
        c = '0';
      }
    }
    Instant const epoch = parse_time(epoch_string);
    read_line();
    for (int i = 0; i < orbits_.size(); ++i) {
      // P record: the position and clock record.
      CHECK_EQ(column(1), 'P') << location;
      SatelliteIdentifier id;
      id.group = version_ == Version::A ? SatelliteGroup::GPS
                                        : SatelliteGroup{column(2)};
      id.index = integer_columns(3, 4);
      auto const it = orbits_.find(id);
      CHECK(it != orbits_.end()) << "Unknown satellite identifier "
                                 << id << ": " << location;

      // The SP3-c and SP3-d specification require that the satellite order of
      // the P, EP, V, and EV records be the same as the order of the satellite
      // ID records.
      // This wording was added to the SP3-c specification by the 2006-09-27
      // amendment, which describes it as a “clarification”, so the intent seems
      // to be that this was required from the start for SP3-c, and perhaps for
      // earlier versions as well.
      // If this breaks for SP3-a or SP3-b, consider exempting these versions
      // from the check.
      CHECK_EQ(id, satellites_[i]) << location;

      std::vector<not_null<std::unique_ptr<DiscreteTrajectory<ITRS>>>>& orbit =
          it->second;
      DiscreteTrajectory<ITRS>& arc = *orbit.back();

      Position<ITRS> const position =
          Displacement<ITRS>({float_columns(5, 18) * Kilo(Metre),
                              float_columns(19, 32) * Kilo(Metre),
                              float_columns(33, 46) * Kilo(Metre)}) +
          ITRS::origin;
      // If the file does not provide velocities, fill the trajectory with NaN
      // velocities; we then replace it with another trajectory whose velocities
      // are computed using a finite difference formula.
      Velocity<ITRS> velocity({NaN<Speed>, NaN<Speed>, NaN<Speed>});

      read_line();
      if (version_ >= Version::C && line.has_value() && columns(1, 2) == "EP") {
        // Ignore the optional EP record (the position and clock correlation
        // record).
        read_line();
      }

      if (has_velocities_) {
        // V record: the velocity and clock rate-of-change record.
        CHECK_EQ(column(1), 'V') << location;
        if (version_ > Version::A) {
          CHECK_EQ(SatelliteGroup{column(2)}, id.group) << location;
        }
        CHECK_EQ(integer_columns(3, 4), id.index) << location;
        Speed const speed_unit =
            dialect == Dialect::GRGS ? Metre / Second : Deci(Metre) / Second;
        velocity = Velocity<ITRS>({float_columns(5, 18) * speed_unit,
                                   float_columns(19, 32) * speed_unit,
                                   float_columns(33, 46) * speed_unit});

        read_line();
        if (version_ >= Version::C && line.has_value() &&
            columns(1, 2) == "EV") {
          // Ignore the optional EV record (the velocity and clock
          // rate-of-change correlation record).
          read_line();
        }
      }

      // Bad or absent positional and velocity values are to be set to 0.000000.
      if (position == ITRS::origin || velocity == ITRS::unmoving) {
        if (!arc.Empty()) {
          orbit.push_back(make_not_null_unique<DiscreteTrajectory<ITRS>>());
        }
      } else {
        arc.Append(epoch, {position, velocity});
      }
    }
  }
  if (dialect != Dialect::ILRSA) {
    CHECK_EQ(columns(1, 3), "EOF") << location;
    read_line();
  }
  for (auto& [id, orbit] : orbits_) {
    // Do not leave a final empty trajectory if the orbit ends with missing
    // data.
    if (orbit.back()->Empty()) {
      orbit.pop_back();
    }
  }
  CHECK(!line.has_value()) << location;
  if (!has_velocities_) {
    for (auto& [id, orbit] : orbits_) {
      for (auto& arc : orbit) {
#define COMPUTE_VELOCITIES_CASE(n)            \
          case n:                             \
            arc = ComputeVelocities<n>(*arc); \
            break

        switch (arc->Size()) {
          COMPUTE_VELOCITIES_CASE(1);
          COMPUTE_VELOCITIES_CASE(2);
          COMPUTE_VELOCITIES_CASE(3);
          COMPUTE_VELOCITIES_CASE(4);
          COMPUTE_VELOCITIES_CASE(5);
          COMPUTE_VELOCITIES_CASE(6);
          COMPUTE_VELOCITIES_CASE(7);
          COMPUTE_VELOCITIES_CASE(8);
          default:
            arc = ComputeVelocities<9>(*arc);
            break;
        }

#undef COMPUTE_VELOCITIES_CASE
      }
    }
  }
  for (auto const& [id, orbit] : orbits_) {
    auto const [it, inserted] =
        const_orbits_.emplace(std::piecewise_construct,
                              std::forward_as_tuple(id),
                              std::forward_as_tuple());
    CHECK(inserted) << id;
    auto& const_orbit = it->second;
    for (auto const& arc : orbit) {
      const_orbit.push_back(arc.get());
    }
  }
}

std::vector<StandardProduct3::SatelliteIdentifier> const&
StandardProduct3::satellites() const {
  return satellites_;
}

std::vector<not_null<DiscreteTrajectory<ITRS> const*>> const&
StandardProduct3::orbit(SatelliteIdentifier const& id) const {
  return FindOrDie(const_orbits_, id);
}

StandardProduct3::Version StandardProduct3::version() const {
  return version_;
}

bool StandardProduct3::file_has_velocities() const {
  return has_velocities_;
}

bool operator==(StandardProduct3::SatelliteIdentifier const& left,
                StandardProduct3::SatelliteIdentifier const& right) {
  return left.group == right.group && left.index == right.index;
}

bool operator<(StandardProduct3::SatelliteIdentifier const& left,
               StandardProduct3::SatelliteIdentifier const& right) {
  return left.group < right.group ||
         (left.group == right.group && left.index < right.index);
}

std::ostream& operator<<(std::ostream& out,
                         StandardProduct3::Version const& version) {
  return out << std::string(1, static_cast<char>(version));
}

std::ostream& operator<<(std::ostream& out,
                         StandardProduct3::Dialect const& dialect) {
  switch (dialect) {
    case StandardProduct3::Dialect::Standard:
      return out << "Standard SP3";
    case StandardProduct3::Dialect::ILRSA:
      return out << "ILRSA SP3 dialect";
    case StandardProduct3::Dialect::ILRSB:
      return out << "ILRSB SP3 dialect";
    case StandardProduct3::Dialect::GRGS:
      return out << "GRGS SP3 dialect";
    case StandardProduct3::Dialect::ChineseMGEX:
      return out << "SHAO and WHU MGEX SP3 dialect";
    default:
      return out << "Unknown SP3 dialect (" << static_cast<int>(dialect) << ")";
  }
}

std::ostream& operator<<(std::ostream& out,
                         StandardProduct3::SatelliteGroup const& group) {
  return out << std::string(1, static_cast<char>(group));
}

std::ostream& operator<<(std::ostream& out,
                         StandardProduct3::SatelliteIdentifier const& id) {
  return out << id.group << absl::StrCat(absl::Dec(id.index, absl::kZeroPad2));
}

}  // namespace internal_standard_product_3
}  // namespace astronomy
}  // namespace principia
