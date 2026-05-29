#include "ksp_plugin/plugin_reader.hpp"

#include <memory>
#include <ranges>
#include <span>
#include <string>
#include <utility>

#include "absl/log/log_sink_registry.h"
#include "base/graveyard.hpp"
#include "base/stoppable_thread.hpp"

namespace principia {
namespace ksp_plugin {
namespace _plugin_reader {
namespace internal {

using namespace principia::base::_graveyard;
using namespace principia::base::_stoppable_thread;

// This is used to destroy a single arena after loading the plugin, so there is
// no need for multiple threads.
auto* const gravedigger = new Graveyard(1);

PluginReader::PluginReader(
    serialization::Plugin const& message,
    not_null<std::unique_ptr<google::protobuf::Arena>> arena)
    : reader_(MakeStoppableThread(
          [this, arena = std::move(arena), &message]() mutable {
            absl::AddLogSink(&log_sink_);
            auto result =
                Plugin::ReadFromMessage(message, [this](bool will_be_slow) {
                  absl::MutexLock l(lock_);
                  will_be_slow_ = will_be_slow;
                });
            gravedigger->Bury(std::move(arena));
            absl::RemoveLogSink(&log_sink_);
            absl::MutexLock l(lock_);
            will_be_slow_ = false;
            result_ = std::move(result);
          })) {}

bool PluginReader::WillBeSlow() const {
  auto slowness_known = [this]() {
    lock_.AssertReaderHeld();
    return will_be_slow_.has_value();
  };
  absl::ReaderMutexLock l(lock_);
  lock_.Await(absl::Condition(&slowness_known));
  return *will_be_slow_;
}

not_null<std::unique_ptr<Plugin>> PluginReader::Await() {
  auto has_result = [this]() {
    lock_.AssertReaderHeld();
    return result_ != nullptr;
  };
  absl::MutexLock l(lock_);
  lock_.Await(absl::Condition(&has_result));
  return std::move(result_);
}

std::unique_ptr<Plugin> PluginReader::get() {
  absl::MutexLock l(lock_);
  return std::move(result_);
}

std::string const& PluginReader::logs() {
  auto const lines = log_sink_.logs();
  std::span tail = lines;
  if (tail.size() > 10) {
    tail = tail.subspan(tail.size() - 10);
  }
  logs_snapshot_ =
      tail | std::ranges::views::join | std::ranges::to<std::string>();
  return logs_snapshot_;
}

}  // namespace internal
}  // namespace _plugin_reader
}  // namespace ksp_plugin
}  // namespace principia
