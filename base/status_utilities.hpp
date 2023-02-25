#pragma once

#include "absl/status/status.h"
#include "absl/status/statusor.h"

namespace principia {
namespace base {
namespace _status_utilities {
namespace internal {

inline absl::Status const& GetStatus(absl::Status const& s) {
  return s;
}

template<typename T>
absl::Status const& GetStatus(absl::StatusOr<T> const& s) {
  return s.status();
}

}  // namespace internal

#define CHECK_OK(value) CHECK_EQ((value), ::absl::OkStatus())
#define DCHECK_OK(value) DCHECK_EQ((value), ::absl::OkStatus())

#define RETURN_IF_ERROR(expr)                                                \
  do {                                                                       \
    /* Using _status below to avoid capture problems if expr is "status". */ \
    ::absl::Status const _status =                                           \
        (::principia::base::internal_status_utilities::GetStatus(expr));     \
    if (!_status.ok())                                                       \
      return _status;                                                        \
  } while (false)

}  // namespace _status_utilities
}  // namespace base
}  // namespace principia

namespace principia::base {
using namespace principia::base::_status_utilities;
}  // namespace principia::base
