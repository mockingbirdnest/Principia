#pragma once

#include "mathematica/logger.hpp"

#include <string>
#include <utility>

#include "mathematica/mathematica.hpp"

namespace principia {
namespace mathematica {
namespace internal_logger {

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
  for (auto const& [name, values] : name_and_multiple_values_) {
    file_ << internal_mathematica::RawApply(
                 "Set",
                 {name, internal_mathematica::RawApply("List", values)}) +
                 ";\n";
  }
  for (auto const& [name, value] : name_and_single_value_) {
    file_ << internal_mathematica::RawApply("Set", {name, value}) + ";\n";
  }
}

template<typename... Args>
void Logger::Append(std::string const& name, Args... args) {
  if (enabled_) {
    name_and_multiple_values_[name].push_back(ToMathematica(args...));
  }
}

template<typename... Args>
void Logger::Set(std::string const& name, Args... args) {
  if (enabled_) {
    name_and_single_value_[name] = ToMathematica(args...);
  }
}

inline void Logger::Enable() {
  enabled_ = true;
}

inline void Logger::Disable() {
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

inline std::atomic_uint64_t Logger::id_ = 0;
inline Logger::ConstructionCallback Logger::construction_callback_ = nullptr;
inline ABSL_CONST_INIT absl::Mutex Logger::construction_callback_lock_(
    absl::kConstInit);

}  // namespace internal_logger
}  // namespace mathematica
}  // namespace principia
