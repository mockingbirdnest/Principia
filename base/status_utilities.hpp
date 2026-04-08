#pragma once

#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "base/status_utilities.hpp"  // 🧙 For RETURN_IF_ERROR.

namespace principia {
namespace base {
namespace _status_utilities {
namespace internal {

inline absl::Status GetStatus(absl::Status const& s) {
  return s;
}

template<typename T>
absl::Status GetStatus(absl::StatusOr<T> const& s) {
  return s.status();
}

}  // namespace internal

#define RETURN_IF_ERROR(expr)                                                \
  do {                                                                       \
    /* Using _status below to avoid capture problems if expr is "status". */ \
    /* Not const to allow moving. */                                         \
    ::absl::Status _status =                                                 \
        (::principia::base::_status_utilities::internal::GetStatus(expr));   \
    if (!_status.ok())                                                       \
      return _status;                                                        \
  } while (false)

}  // namespace _status_utilities
}  // namespace base
}  // namespace principia
