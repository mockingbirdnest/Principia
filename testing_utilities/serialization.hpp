#pragma once

#include <cstdint>
#include <filesystem>
#include <string>
#include <vector>

#include "base/array.hpp"
#include "base/macros.hpp"
#include "serialization/testing_utilities.pb.h"

namespace principia {
namespace testing_utilities {
namespace _serialization {
namespace internal {

using namespace principia::base::_array;

#if PRINCIPIA_COMPILER_MSVC
inline std::u16string ReadFromBase32768File(
    std::filesystem::path const& filename);
#endif

inline std::vector<std::uint8_t> ReadFromBinaryFile(
    std::filesystem::path const& filename);

inline std::string ReadFromHexadecimalFile(
    std::filesystem::path const& filename);

inline std::vector<std::string> ReadLinesFromBase64File(
    std::filesystem::path const& filename);

inline std::vector<std::string> ReadLinesFromHexadecimalFile(
    std::filesystem::path const& filename);

inline serialization::TabulatedData ReadFromTabulatedData(
    std::filesystem::path const& filename);

#if PRINCIPIA_COMPILER_MSVC
inline void WriteToBase32768File(std::filesystem::path const& filename,
                                 Array<std::uint8_t const> serialized);
#endif

inline void WriteToBinaryFile(std::filesystem::path const& filename,
                              Array<std::uint8_t const> serialized);

inline void WriteToHexadecimalFile(std::filesystem::path const& filename,
                                   Array<std::uint8_t const> serialized);

}  // namespace internal

#if PRINCIPIA_COMPILER_MSVC
using internal::ReadFromBase32768File;
#endif
using internal::ReadFromBinaryFile;
using internal::ReadFromHexadecimalFile;
using internal::ReadFromTabulatedData;
using internal::ReadLinesFromBase64File;
using internal::ReadLinesFromHexadecimalFile;
#if PRINCIPIA_COMPILER_MSVC
using internal::WriteToBase32768File;
#endif
using internal::WriteToBinaryFile;
using internal::WriteToHexadecimalFile;

}  // namespace _serialization
}  // namespace testing_utilities
}  // namespace principia

namespace principia::testing_utilities {
using namespace principia::testing_utilities::_serialization;
}  // namespace principia::testing_utilities

#include "testing_utilities/serialization_body.hpp"
