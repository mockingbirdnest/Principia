#pragma once

#include <experimental/filesystem>
#include <string>

namespace principia {
namespace testing_utilities {

std::string ReadFromBinaryFile(
    std::experimental::filesystem::path const& filename);

std::string ReadFromHexadecimalFile(
    std::experimental::filesystem::path const& filename);

void WriteToBinaryFile(std::experimental::filesystem::path const& filename,
                       std::string const& serialized);

void WriteToHexadecimalFile(std::experimental::filesystem::path const& filename,
                            std::string const& serialized);

}  // namespace testing_utilities
}  // namespace principia

#include "testing_utilities/serialization_body.hpp"
