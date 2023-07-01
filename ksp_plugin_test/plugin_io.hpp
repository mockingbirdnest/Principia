#pragma once

#include <filesystem>
#include <memory>
#include <string_view>

#include "base/not_null.hpp"
#include "ksp_plugin/plugin.hpp"

namespace principia {
namespace ksp_plugin_test {
namespace _plugin_io {
namespace internal {

using namespace principia::base::_not_null;
using namespace principia::ksp_plugin::_plugin;

// Reads a plugin from a file containing only the "serialized_plugin = " lines,
// with "serialized_plugin = " dropped.
not_null<std::unique_ptr<Plugin const>> ReadPluginFromFile(
    std::filesystem::path const& filename,
    std::string_view compressor,
    std::string_view encoder);

// Same as above, but records the number of bytes read from the file.
not_null<std::unique_ptr<Plugin const>> ReadPluginFromFile(
    std::filesystem::path const& filename,
    std::string_view compressor,
    std::string_view encoder,
    std::int64_t& bytes_processed);

// Writes a plugin to a file.
void WritePluginToFile(std::filesystem::path const& filename,
                       std::string_view compressor,
                       std::string_view encoder,
                       not_null<std::unique_ptr<Plugin const>> plugin);

// Same as above, but records the number of bytes written to the file.
void WritePluginToFile(std::filesystem::path const& filename,
                       std::string_view compressor,
                       std::string_view encoder,
                       not_null<std::unique_ptr<Plugin const>> plugin,
                       std::int64_t& bytes_processed);

}  // namespace internal

using internal::ReadPluginFromFile;
using internal::WritePluginToFile;

}  // namespace _plugin_io
}  // namespace ksp_plugin_test
}  // namespace principia
