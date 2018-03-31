
#pragma once

#include <optional>

// We need to pollute the |::| namespace in order for lookup to work, like
// glog/stl_logging.h.  The following caveat from glog/stl_logging.h applies
// here too.

// Note that if you want to use these operators from the non-global namespace,
// you may get an error since they are not in namespace std (and they are not
// in namespace std since that would result in undefined behavior).  You may
// need to write
//
//   using ::operator<<;
//
// to fix these errors.

// If |optional|, logs |*optional|, otherwise, logs |"nullopt"|.
template<typename T>
std::ostream& operator<<(std::ostream& out,
                         std::optional<T> const& optional);

#include "base/optional_logging_body.hpp"
