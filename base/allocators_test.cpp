
#include "base/allocators.hpp"

#include <vector>

#include "gmock/gmock.h"
#include "gtest/gtest.h"

using testing::Each;

namespace principia {
namespace base {

namespace {
template <typename T>
using default_vector = std::vector<T, DefaultInitializationAllocator<T>>;
}  //namespace

class AllocatorsTest : public testing::Test {};

TEST_F(AllocatorsTest, Resize) {
  default_vector<uint8_t> default_initialized;
  std::size_t const size = 10000;
  default_initialized.reserve(size);
  // Possibly UB, but for practical purposes it should work (and this is a test
  // anyway).
  memset(default_initialized.data(), 42, size);
  default_initialized.resize(size);
  EXPECT_THAT(default_initialized, Each('*'));
  std::vector<uint8_t> value_initialized;
  value_initialized.reserve(size);
  memset(value_initialized.data(), 42, size);
  value_initialized.resize(size);
  EXPECT_THAT(value_initialized, Each('\0'));
}

// TODO(egg): is there a way to test the constructor that doesn't involve access
// violations?

}  // namespace base
}  // namespace principia
