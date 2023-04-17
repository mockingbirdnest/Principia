#pragma once

#include <atomic>
#include <filesystem>
#include <functional>
#include <map>
#include <string>
#include <vector>

#include "absl/synchronization/mutex.h"
#include "base/file.hpp"
#include "base/not_null.hpp"

namespace principia {
namespace mathematica {
namespace _logger {
namespace internal {

using namespace principia::base::_file;
using namespace principia::base::_not_null;

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

  // When a logger is disabled, the calls to |Append| and |Set| have no effect.
  // Loggers are enabled at construction.
  void Enable();
  void Disable();

  // It is possible to register a (single) callback that gets executed at
  // construction of each logger.  This makes it possible for distant code to
  // perform operations on the logger, e.g., to disable it or log client-
  // specific data.  The path passed to the callback is the path passed to the
  // constructor (i.e., before any alteration performed by |make_unique|).  The
  // |id| is the uniquely-generated id for the logger (if |make_unique| is true)
  // or nullopt (if |make_unique| is false).
  using ConstructionCallback =
      std::function<void(std::filesystem::path const&,
                         std::optional<std::uint64_t> id,
                         not_null<Logger*>)>;

  static void SetConstructionCallback(ConstructionCallback callback);
  static void ClearConstructionCallback();

 private:
  std::atomic_bool enabled_ = true;
  std::optional<std::uint64_t> my_id_;
  OFStream file_;
  std::map<std::string, std::vector<std::string>> name_and_multiple_values_;
  std::map<std::string, std::string> name_and_single_value_;

  static std::atomic_uint64_t id_;
  static ConstructionCallback construction_callback_;
  static absl::Mutex construction_callback_lock_;
};

}  // namespace internal

using internal::Logger;

}  // namespace _logger
}  // namespace mathematica
}  // namespace principia

// It is not strictly necessary to have the body inline, but it avoids having to
// add the .cpp file to a project when using the logger.
#include "mathematica/logger_body.hpp"
