#pragma once

#include "mathematica/logger.hpp"

#include <string>

#include "mathematica/mathematica.hpp"

namespace principia {
namespace mathematica {
namespace internal_logger {

inline Logger::Logger(std::filesystem::path const& path, bool const make_unique)
    : file_([make_unique, &path]() {
        if (make_unique || PRINCIPIA_MATHEMATICA_LOGGER_REGRESSION_TEST != 0) {
          std::filesystem::path filename = path.stem();
          if (make_unique) {
            filename += std::to_string(id_++);
          }
#if PRINCIPIA_MATHEMATICA_LOGGER_REGRESSION_TEST
          filename += "_new";
#endif
          filename += path.extension();
          return path.parent_path() / filename;
        } else {
          return path;
        }
      }()) {}

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
  name_and_multiple_values_[name].push_back(ToMathematica(args...));
}

template<typename... Args>
void Logger::Set(std::string const& name, Args... args) {
  name_and_single_value_[name] = ToMathematica(args...);
}

inline std::atomic_uint64_t Logger::id_ = 0;

}  // namespace internal_logger
}  // namespace mathematica
}  // namespace principia
