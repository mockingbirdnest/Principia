
#include "base/allocator_deleter.hpp"

#include <memory>

#include "base/allocated_by.hpp"
#include "base/allocator_new.hpp"
#include "base/malloc_allocator.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

namespace principia {
namespace base {

using testing::IsNull;
using testing::Not;
using testing::Pointee;

// Class which counts how many times its constructor and destructor have been
// called.
class CtorDtorCounter {
 public:
  CtorDtorCounter(int& ctor_count, int& dtor_count) : dtor_count_(&dtor_count) {
    ++ctor_count;
  }

  ~CtorDtorCounter() {
    ++*dtor_count_;
  }

 private:
  int* dtor_count_;
};

TEST(AllocatorDeleterTest, Lifecycle) {
  // Create object.
  int ctor_count = 0;
  int dtor_count = 0;
  CtorDtorCounter* p = AssertAllocatedBy<std::allocator<CtorDtorCounter>>(
      new (AllocateWith<std::allocator<CtorDtorCounter>>{})
          CtorDtorCounter(ctor_count, dtor_count));
  EXPECT_EQ(ctor_count, 1);
  EXPECT_EQ(dtor_count, 0);

  AllocatorDeleter<std::allocator<CtorDtorCounter>> deleter;
  deleter(AssertAllocatedBy<std::allocator<CtorDtorCounter>>(p));
  EXPECT_EQ(ctor_count, 1);
  EXPECT_EQ(dtor_count, 1);
}

TEST(AllocatorDeleterTest, MemberTypes) {
  EXPECT_TRUE((std::is_same<AllocatorDeleter<std::allocator<int>>::pointer,
                            AllocatedBy<std::allocator<int>>>::value));
  EXPECT_TRUE((std::is_same<AllocatorDeleter<std::allocator<int>>::value_type,
                            int>::value));
}

TEST(AllocatorDeleterTest, Upcasting) {
  class Base {};
  class Sub : public Base {};

  AllocatorDeleter<std::allocator<Sub>> a;

  AllocatorDeleter<std::allocator<Base>> b(a);
}

TEST(AllocatorDeleterTest, MakeUniqueWithAllocator) {
  std::unique_ptr<int, AllocatorDeleter<std::allocator<int>>> owned =
      MakeUniqueWithAllocator<std::allocator<int>>(2);

  EXPECT_THAT(owned, Pointee(2));
}

TEST(AllocatorDeleterTest, ShareUnique) {
  auto owned = MakeUniqueWithAllocator<std::allocator<int>>(2);
  std::shared_ptr<int> shared = std::move(owned);

  EXPECT_THAT(shared, Pointee(2));
}

}  // namespace base
}  // namespace principia
