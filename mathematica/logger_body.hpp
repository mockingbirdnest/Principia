#pragma once

#include "mathematica/logger.hpp"

#include <string>
#include <utility>

#include "mathematica/mathematica.hpp"

namespace principia {
namespace mathematica {
namespace _logger {
namespace internal {

using namespace principia::mathematica::_mathematica;

inline Logger::Logger(std::filesystem::path const& path, bool const make_unique)
    : file_([this, make_unique, &path]() {
        if (make_unique || PRINCIPIA_MATHEMATICA_LOGGER_REGRESSION_TEST != 0) {
          std::filesystem::path filename = path.stem();
          if (make_unique) {
            my_id_ = id_++;
            filename += std::to_string(my_id_.value());
          }
#if PRINCIPIA_MATHEMATICA_LOGGER_REGRESSION_TEST
          filename += "_new";
#endif
          filename += path.extension();
          return path.parent_path() / filename;
        } else {
          return path;
        }
      }()) {
  absl::ReaderMutexLock l(&construction_callback_lock_);
  if (construction_callback_ != nullptr) {
    construction_callback_(path, my_id_, this);
  }
}

inline Logger::~Logger() {
  Flush();
}

inline void Logger::Flush() {
  absl::MutexLock l(&lock_);
  FlushLocked();
}

inline void Logger::FlushAndClear() {
  absl::MutexLock l(&lock_);
  FlushLocked();
  name_and_multiple_values_.clear();
  name_and_single_value_.clear();
}

template<typename... Args>
void Logger::Append(std::string const& name, Args... args) {
  absl::MutexLock l(&lock_);
  if (enabled_) {
    name_and_multiple_values_[name].push_back(ToMathematica(args...));
  }
}

template<typename... Args>
void Logger::Set(std::string const& name, Args const&... args) {
  absl::MutexLock l(&lock_);
  if (enabled_) {
    name_and_single_value_[name] = ToMathematica(args...);
  }
}

inline void Logger::Enable() {
  absl::MutexLock l(&lock_);
  enabled_ = true;
}

inline void Logger::Disable() {
  absl::MutexLock l(&lock_);
  enabled_ = false;
}

inline void Logger::SetConstructionCallback(ConstructionCallback callback) {
  absl::MutexLock l(&construction_callback_lock_);
  construction_callback_ = std::move(callback);
}

inline void Logger::ClearConstructionCallback() {
  absl::MutexLock l(&construction_callback_lock_);
  construction_callback_ = nullptr;
}

inline void Logger::FlushLocked() {
  lock_.AssertHeld();
  using _mathematica::internal::RawApply;
  for (auto const& [name, values] : name_and_multiple_values_) {
    file_ << RawApply("Set", {name, RawApply("List", values)}) + ";\n";
  }
  for (auto const& [name, value] : name_and_single_value_) {
    file_ << RawApply("Set", {name, value}) + ";\n";
  }
  file_.Flush();
}

inline std::atomic_uint64_t Logger::id_ = 0;
inline Logger::ConstructionCallback Logger::construction_callback_ = nullptr;
ABSL_CONST_INIT inline absl::Mutex Logger::construction_callback_lock_(
    absl::kConstInit);

}  // namespace internal
}  // namespace _logger
}  // namespace mathematica
}  // namespace principia
