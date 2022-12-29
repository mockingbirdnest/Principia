#pragma once

#include <atomic>
#include <filesystem>
#include <map>
#include <string>
#include <vector>

#include "base/file.hpp"

namespace principia {
namespace mathematica {
namespace internal_logger {

using base::OFStream;

// Define this value to 1 to force the logger to append "_new" to the file
// names, which is useful for regression testing of the logger.
#define PRINCIPIA_MATHEMATICA_LOGGER_REGRESSION_TEST 0

// An RAII object to help with Mathematica logging.
class Logger final {
 public:
  // Creates a logger object that will, at destruction, write to the given file.
  // If make_unique is true, a unique id is inserted before the file extension
  // to identify different loggers.
  Logger(std::filesystem::path const& path,
         bool make_unique = true);
  ~Logger();

  // Flushes the contents of the logger to the file.  This should not normally
  // be called explicitly, but it's useful when debugging crashes.  Beware, the
  // resulting file may contain many assignments to the same variable.
  void Flush();

  // Appends an element to the list of values for the List variable |name|.  The
  // |args...| are passed verbatim to ToMathematica for stringification.  When
  // this object is destroyed, an assignment is generated for each of the
  // variables named in a call to Append.
  template<typename... Args>
  void Append(std::string const& name, Args... args);

  // Sets an element as the single value for the variable |name|.  The
  // |args...| are passed verbatim to ToMathematica for stringification.  When
  // this object is destroyed, an assignment is generated for each of the
  // variables named in a call to Set.
  template<typename... Args>
  void Set(std::string const& name, Args... args);

 private:
  OFStream file_;
  std::map<std::string, std::vector<std::string>> name_and_multiple_values_;
  std::map<std::string, std::string> name_and_single_value_;

  static std::atomic_uint64_t id_;
};

}  // namespace internal_logger

using internal_logger::Logger;

}  // namespace mathematica
}  // namespace principia

// It is not strictly necessary to have the body inline, but it avoids having to
// add the .cpp file to a project when using the logger.
#include "mathematica/logger_body.hpp"
