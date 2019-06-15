#pragma once

#include "testing_utilities/serialization.hpp"

#include <cstdint>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>

#include "base/base32768.hpp"
#include "google/protobuf/io/zero_copy_stream_impl.h"
#include "google/protobuf/text_format.h"
#include "serialization/testing_utilities.pb.h"

namespace principia {
namespace testing_utilities {

inline void PrintUtf16ToFile(char16_t const c, std::fstream& file) {
  char chars[sizeof(char16_t)];
  std::memcpy(chars, &c, sizeof(char16_t));
  file << chars[0] << chars[1];
}

#if PRINCIPIA_COMPILER_MSVC
std::u16string ReadFromBase32768File(
    std::filesystem::path const& filename) {
  std::vector<std::uint8_t> const binary = ReadFromBinaryFile(filename);
  std::u16string base32768(binary.size() >> 1, u'\0');
  std::memcpy(base32768.data(), binary.data(), binary.size());

  // Base 32768 doesn't use ASCII characters so strip them.
  std::u16string result;
  for (char16_t c : base32768) {
    if (c > u'\x00FF' && c != u'\xFEFF') {
      result.append(1, c);
    }
  }
  return result;
}
#endif

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

std::vector<std::string> ReadLinesFromHexadecimalFile(
    std::filesystem::path const& filename) {
  std::fstream file = std::fstream(filename);
  CHECK(file.good()) << filename;
  std::vector<std::string> hex;
  while (!file.eof()) {
    std::string line;
    std::getline(file, line);
    hex.push_back("");
    for (auto const c : line) {
      if ((c >= '0' && c <= '9') || (c >= 'A' && c <= 'F')) {
        hex.back().push_back(c);
      }
    }
    if (hex.back().empty()) {
      hex.pop_back();
    }
  }
  file.close();
  return hex;
}

serialization::TabulatedData ReadFromTabulatedData(
    std::filesystem::path const& filename) {
  serialization::TabulatedData tabulated_data;
  std::ifstream tabulated_data_ifstream(filename);
  CHECK(tabulated_data_ifstream.good());
  google::protobuf::io::IstreamInputStream tabulated_data_zcs(
      &tabulated_data_ifstream);
  CHECK(google::protobuf::TextFormat::Parse(&tabulated_data_zcs,
                                            &tabulated_data));
  CHECK_LT(0, tabulated_data.entry_size());
  return tabulated_data;
}

#if PRINCIPIA_COMPILER_MSVC
void WriteToBase32768File(std::filesystem::path const& filename,
                          base::Array<std::uint8_t const> serialized) {
  std::fstream file = std::fstream(filename,
                                   std::ios::binary | std::ios::out);
  CHECK(file.good()) << filename;
  base::Base32768Encoder</*null_terminated=*/false> encoder;
  auto const base32768 = encoder.Encode(serialized);

  PrintUtf16ToFile(u'\uFEFF', file);
  int index = 0;
  while (index < base32768.size) {
    PrintUtf16ToFile(base32768.data[index], file);
    ++index;
    if (index % 40 == 0) {
      PrintUtf16ToFile(u'\r', file);
      PrintUtf16ToFile(u'\n', file);
    }
  }
  if (index % 40 != 0) {
      PrintUtf16ToFile(u'\r', file);
      PrintUtf16ToFile(u'\n', file);
  }
  file.close();
}
#endif

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
  int index = 0;
  while (index < serialized.size) {
    file << std::hex << std::uppercase << std::setw(2) << std::setfill('0')
         << static_cast<int>(serialized.data[index]);
    ++index;
    if (index % 40 == 0) {
      file << '\n';
    }
  }
  if (index % 40 != 0) {
    file << '\n';
  }
  file.close();
}

}  // namespace testing_utilities
}  // namespace principia

