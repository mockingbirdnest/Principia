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
                                 base::Array<std::uint8_t const> serialized);
#endif

inline void WriteToBinaryFile(std::filesystem::path const& filename,
                              base::Array<std::uint8_t const> serialized);

inline void WriteToHexadecimalFile(std::filesystem::path const& filename,
                                   base::Array<std::uint8_t const> serialized);

}  // namespace testing_utilities
}  // namespace principia

#include "testing_utilities/serialization_body.hpp"
