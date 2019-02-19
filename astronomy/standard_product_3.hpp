#pragma once

#include <filesystem>
#include <string>

namespace principia {
namespace astronomy {

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
};

}  // namespace astronomy
}  // namespace principia
