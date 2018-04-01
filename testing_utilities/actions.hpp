#pragma once

#include "gmock/gmock.h"

namespace principia {
namespace testing_utilities {

ACTION_TEMPLATE(FillUniquePtr,
                // Note the comma between int and k:
                HAS_1_TEMPLATE_PARAMS(int, k),
                AND_1_VALUE_PARAMS(ptr)) {
  std::get<k>(args)->reset(ptr);
}

}  // namespace testing_utilities
}  // namespace principia
