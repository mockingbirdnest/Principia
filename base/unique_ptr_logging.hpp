#pragma once

#include <memory>

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

// Logs the pointer.  No transfer of ownership.
template<typename T, typename U>
std::ostream& operator<<(std::ostream& out,
                         std::unique_ptr<T, U> const& pointer);

#include "base/unique_ptr_logging_body.hpp"
