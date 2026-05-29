#pragma once

#include <future>
#include <functional>
#include <memory>
#include <string>

#include "absl/synchronization/mutex.h"
#include "base/not_null.hpp"
#include "base/string_log_sink.hpp"
#include "google/protobuf/arena.h"
#include "ksp_plugin/plugin.hpp"

namespace principia {
namespace ksp_plugin {
namespace _plugin_reader {
namespace internal {

using namespace principia::base::_not_null;
using namespace principia::base::_string_log_sink;
using namespace principia::ksp_plugin::_plugin;

class PluginReader {
 public:
  // `message` must be allocated on `*arena`.  The arena will be destroyed at
  // some point after this constructor is called, possibly after the destruction
  // of the constructed object (the arena is sent to a graveyard).
  PluginReader(
      serialization::Plugin const& message,
      not_null<std::unique_ptr<google::protobuf::Arena>> arena);

  // Blocks until `Plugin::ReadFromMessage` has determined whether
  // deserialization will be slow, and returns the result of that determination.
  bool WillBeSlow() const;

  // Waits for `Plugin::ReadFromMessage` to finish and returns its result.
  not_null<std::unique_ptr<Plugin>> Await();
  // Returns the result of `Plugin::ReadFromMessage` if available, and `nullptr`
  // otherwise.
  std::unique_ptr<Plugin> get();

  // Returns the last few log lines.
  // The result is owned by this object and is is not mutated until the next
  // call to `logs()`.
  std::string const& logs();

 private:
  StringLogSink log_sink_;
  std::string logs_snapshot_;
  mutable absl::Mutex lock_;
  std::optional<bool> will_be_slow_ ABSL_GUARDED_BY(lock_);
  std::unique_ptr<Plugin> result_ ABSL_GUARDED_BY(lock_);

  std::jthread reader_;
};

}  // namespace internal

using internal::PluginReader;

}  // namespace _plugin_reader
}  // namespace ksp_plugin
}  // namespace principia
