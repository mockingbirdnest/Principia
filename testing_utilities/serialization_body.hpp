#pragma once

#include "testing_utilities/serialization.hpp"

#include <fstream>
#include <iomanip>
#include <string>

namespace principia {
namespace testing_utilities {

// TODO(egg): This should use |UniqueBytes| or |Bytes| rather than |std::string|
// for non-text bytes.

inline std::string ReadFromBinaryFile(std::filesystem::path const& filename) {
  std::fstream file = std::fstream(filename, std::ios::in | std::ios::binary);
  CHECK(file.good()) << filename;
  std::string binary;
  while (!file.eof()) {
    char c;
    file.get(c);
    binary.append(1, c);
  }
  file.close();
  return binary;
}

inline std::string ReadFromHexadecimalFile(
    std::filesystem::path const& filename) {
  std::fstream file = std::fstream(filename);
  CHECK(file.good()) << filename;
  std::string hex;
  while (!file.eof()) {
    std::string line;
    std::getline(file, line);
    for (auto const c : line) {
      if ((c >= '0' && c <= '9') || (c >= 'A' && c <= 'F')) {
        hex.append(1, c);
      }
    }
  }
  file.close();
  return hex;
}

inline void WriteToBinaryFile(std::filesystem::path const& filename,
                              std::string const& serialized) {
  std::fstream file = std::fstream(filename,
                                   std::ios::binary | std::ios::out);
  CHECK(file.good()) << filename;
  file.write(serialized.c_str(), serialized.size());
  file.close();
}

inline void WriteToHexadecimalFile(std::filesystem::path const& filename,
                                   std::string const& serialized) {
  std::fstream file = std::fstream(filename,
                                   std::ios::out);
  CHECK(file.good()) << filename;
  int index = 0;
  for (unsigned char c : serialized) {
    file << std::hex << std::uppercase << std::setw(2) << std::setfill('0')
         << static_cast<int>(c);
    ++index;
    if (index == 40) {
      file << '\n';
      index = 0;
    }
  }
  if (index != 0) {
    file << '\n';
  }
  file.close();
}

}  // namespace testing_utilities
}  // namespace principia

