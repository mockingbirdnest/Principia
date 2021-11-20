#pragma once

#include "absl/status/status.h"
#include "absl/status/statusor.h"

namespace principia {
namespace base {
namespace internal_status_utilities {

inline absl::Status const& GetStatus(absl::Status const& s) {
  return s;
}

template<typename T>
absl::Status const& GetStatus(absl::StatusOr<T> const& s) {
  return s.status();
}

}  // namespace internal_status_utilities

#define CHECK_OK(value) CHECK_EQ((value), ::absl::OkStatus())
#define DCHECK_OK(value) DCHECK_EQ((value), ::absl::OkStatus())
#define EXPECT_OK(value) \
  EXPECT_THAT((value), ::principia::testing_utilities::IsOk());

#define RETURN_IF_ERROR(expr)                                                \
  do {                                                                       \
    /* Using _status below to avoid capture problems if expr is "status". */ \
    ::absl::Status const _status =                                           \
        (::principia::base::internal_status_utilities::GetStatus(expr));     \
    if (!_status.ok())                                                       \
      return _status;                                                        \
  } while (false)

}  // namespace base
}  // namespace principia
