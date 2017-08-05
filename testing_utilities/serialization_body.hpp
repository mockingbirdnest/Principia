#pragma once

#include "testing_utilities/serialization.hpp"

#include <fstream>
#include <iomanip>

namespace principia {
namespace testing_utilities {

std::string ReadFromBinaryFile(
    std::experimental::filesystem::path const& filename) {
  std::fstream file = std::fstream(filename, std::ios::in | std::ios::binary);
  CHECK(file.good());
  std::string binary;
  while (!file.eof()) {
    char c;
    file.get(c);
    binary.append(1, c);
  }
  file.close();
  return binary;
}

std::string ReadFromHexadecimalFile(
    std::experimental::filesystem::path const& filename) {
  std::fstream file = std::fstream(filename);
  CHECK(file.good());
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

void WriteToBinaryFile(std::experimental::filesystem::path const& filename,
                       std::string const& serialized) {
  std::fstream file = std::fstream(filename,
                                   std::ios::binary | std::ios::out);
  CHECK(file.good());
  file.write(serialized.c_str(), serialized.size());
  file.close();
}

void WriteToHexadecimalFile(std::experimental::filesystem::path const& filename,
                            std::string const& serialized) {
  std::fstream file = std::fstream(filename,
                                   std::ios::out);
  CHECK(file.good());
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

