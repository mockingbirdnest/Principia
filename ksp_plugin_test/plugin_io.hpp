#pragma once

#include <filesystem>
#include <memory>
#include <string_view>

#include "base/not_null.hpp"
#include "ksp_plugin/plugin.hpp"

namespace principia {
namespace interface {
namespace internal_plugin_io {

using base::not_null;
using ksp_plugin::Plugin;

// TODO(phl): Use these functions in benchmark.cpp.

// Reads a plugin from a file containing only the "serialized_plugin = " lines,
// with "serialized_plugin = " dropped.
not_null<std::unique_ptr<Plugin const>> ReadPluginFromFile(
    std::filesystem::path const& filename,
    std::string_view compressor,
    std::string_view encoder);

  // Writes a plugin to a file.
void WritePluginToFile(std::filesystem::path const& filename,
                       std::string_view compressor,
                       std::string_view encoder,
                       not_null<std::unique_ptr<Plugin const>> plugin);

}  // namespace internal_plugin_io

using internal_plugin_io::ReadPluginFromFile;
using internal_plugin_io::WritePluginToFile;

}  // namespace interface
}  // namespace principia
