#pragma once

#include "testing_utilities/serialization.hpp"

#include <cstdint>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>

#include "base/base32768.hpp"

namespace principia {
namespace testing_utilities {

std::u16string ReadFromBase32768File(
    std::filesystem::path const& filename) {
  auto file = std::basic_fstream<char16_t>(filename);
  CHECK(file.good()) << filename;
  std::u16string base32768;
  while (!file.eof()) {
    std::u16string line;
    std::getline(file, line);
    for (auto const c : line) {
      base32768.push_back(c);
    }
  }
  file.close();
  return base32768;
}

std::vector<std::uint8_t> ReadFromBinaryFile(
    std::filesystem::path const& filename) {
  std::fstream file = std::fstream(filename, std::ios::in | std::ios::binary);
  CHECK(file.good()) << filename;
  std::vector<std::uint8_t> binary;
  for (char c; file.get(c);) {
    binary.push_back(static_cast<std::uint8_t>(c));
  }
  file.close();
  return binary;
}

std::string ReadFromHexadecimalFile(
    std::filesystem::path const& filename) {
  std::fstream file = std::fstream(filename);
  CHECK(file.good()) << filename;
  std::string hex;
  while (!file.eof()) {
    std::string line;
    std::getline(file, line);
    for (auto const c : line) {
      if ((c >= '0' && c <= '9') || (c >= 'A' && c <= 'F')) {
        hex.push_back(c);
      }
    }
  }
  file.close();
  return hex;
}

void WriteToBase32768File(std::filesystem::path const& filename,
                          base::Array<std::uint8_t const> serialized) {
  auto file = std::basic_fstream<char16_t>(filename, std::ios::out);
  CHECK(file.good()) << filename;
  auto const base32768 =
      base::Base32768Encode(serialized, /*null_terminated=*/false);
  for (int i = 0; i < base32768.size; ++i) {
    file << base32768.data[i];
    if (i % 40 == 0) {
      file << '\n';
    }
  }
  if (base32768.size % 40 != 0) {
    file << '\n';
  }
  file.close();
}

void WriteToBinaryFile(
    std::filesystem::path const& filename,
    base::Array<std::uint8_t const> const serialized) {
  std::fstream file = std::fstream(filename,
                                   std::ios::binary | std::ios::out);
  CHECK(file.good()) << filename;
  file.write(reinterpret_cast<char const*>(serialized.data), serialized.size);
  file.close();
}

void WriteToHexadecimalFile(
    std::filesystem::path const& filename,
    base::Array<std::uint8_t const> const serialized) {
  std::fstream file = std::fstream(filename, std::ios::out);
  CHECK(file.good()) << filename;
  for (int i = 0; i < serialized.size; ++i) {
    file << std::hex << std::uppercase << std::setw(2) << std::setfill('0')
         << static_cast<int>(serialized.data[i]);
    if (i % 40 == 0) {
      file << '\n';
    }
  }
  if (serialized.size % 40 != 0) {
    file << '\n';
  }
  file.close();
}

}  // namespace testing_utilities
}  // namespace principia

