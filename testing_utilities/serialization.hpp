#pragma once

#include <filesystem>
#include <string>

namespace principia {
namespace testing_utilities {

std::string ReadFromBinaryFile(std::filesystem::path const& filename);

std::string ReadFromHexadecimalFile(std::filesystem::path const& filename);

void WriteToBinaryFile(std::filesystem::path const& filename,
                       std::string const& serialized);

void WriteToHexadecimalFile(std::filesystem::path const& filename,
                            std::string const& serialized);

}  // namespace testing_utilities
}  // namespace principia

#include "testing_utilities/serialization_body.hpp"
