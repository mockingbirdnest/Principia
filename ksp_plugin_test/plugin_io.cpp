#include "ksp_plugin_test/plugin_io.hpp"

#include <string>
#include <vector>

#include "base/file.hpp"
#include "base/pull_serializer.hpp"
#include "base/push_deserializer.hpp"
#include "ksp_plugin/interface.hpp"
#include "testing_utilities/serialization.hpp"

namespace principia {
namespace interface {
namespace internal_plugin_io {

const char preferred_compressor[] = "gipfeli";
const char preferred_encoder[] = "base64";

using testing_utilities::ReadLinesFromBase64File;
using testing_utilities::ReadLinesFromHexadecimalFile;
using namespace principia::base::_file;
using namespace principia::base::_pull_serializer;
using namespace principia::base::_push_deserializer;

not_null<std::unique_ptr<Plugin const>> ReadPluginFromFile(
    std::filesystem::path const& filename,
    std::string_view const compressor,
    std::string_view const encoder) {
  std::int64_t bytes_processed;
  return ReadPluginFromFile(filename, compressor, encoder, bytes_processed);
}

not_null<std::unique_ptr<Plugin const>> ReadPluginFromFile(
    std::filesystem::path const& filename,
    std::string_view const compressor,
    std::string_view const encoder,
    std::int64_t& bytes_processed) {
  Plugin const* plugin = nullptr;

  PushDeserializer* deserializer = nullptr;
  auto const lines =
      encoder == "hexadecimal" ? ReadLinesFromHexadecimalFile(filename)
      : encoder == "base64"    ? ReadLinesFromBase64File(filename)
                               : std::vector<std::string>{};
  CHECK(!lines.empty());

  LOG(ERROR) << "Deserialization starting";
  bytes_processed = 0;
  for (std::string const& line : lines) {
    principia__DeserializePlugin(line.c_str(),
                                 &deserializer,
                                 &plugin,
                                 compressor.data(),
                                 encoder.data());
    bytes_processed += line.size();
  }
  principia__DeserializePlugin("",
                               &deserializer,
                               &plugin,
                               compressor.data(),
                               encoder.data());
  LOG(ERROR) << "Deserialization complete";

  return std::unique_ptr<Plugin const>(plugin);
}

void WritePluginToFile(std::filesystem::path const& filename,
                       std::string_view const compressor,
                       std::string_view const encoder,
                       not_null<std::unique_ptr<Plugin const>> plugin) {
  std::int64_t bytes_processed;
  WritePluginToFile(
      filename, compressor, encoder, std::move(plugin), bytes_processed);
}

void WritePluginToFile(std::filesystem::path const& filename,
                       std::string_view const compressor,
                       std::string_view const encoder,
                       not_null<std::unique_ptr<Plugin const>> plugin,
                       std::int64_t& bytes_processed) {
  OFStream file(filename);
  PullSerializer* serializer = nullptr;
  char const* b64 = nullptr;

  LOG(ERROR) << "Serialization starting";
  bytes_processed = 0;
  for (;;) {
    b64 = principia__SerializePlugin(plugin.get(),
                                     &serializer,
                                     preferred_compressor,
                                     preferred_encoder);
    if (b64 == nullptr) {
      break;
    }
    bytes_processed += std::strlen(b64);
    file << b64 << "\n";
    principia__DeleteString(&b64);
  }
  LOG(ERROR) << "Serialization complete";

  Plugin const* released_plugin = plugin.release();
  principia__DeletePlugin(&released_plugin);
}

}  // namespace internal_plugin_io
}  // namespace interface
}  // namespace principia
