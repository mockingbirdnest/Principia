#include "base/box.hpp"

#include <string>
#include <utility>

#include "gmock/gmock.h"
#include "gtest/gtest.h"

namespace principia {
namespace base {

using ::testing::_;
using ::testing::IsNull;
using ::testing::Not;
using ::testing::NotNull;
using ::testing::Pointee;

// Class which counts how many times its constructor and destructor were called.
class CtorDtorCounter {
 public:
  CtorDtorCounter(int& ctor_count, int& dtor_count) : dtor_count_(&dtor_count) {
    ++ctor_count;
  }

  ~CtorDtorCounter() { ++*dtor_count_; }

 private:
  int* dtor_count_;
};

TEST(BoxTest, Lifecycle) {
  int ctor_count = 0;
  int dtor_count = 0;

  {
    Box<CtorDtorCounter> box(ctor_count, dtor_count);

    EXPECT_EQ(ctor_count, 1);
    EXPECT_EQ(dtor_count, 0);
  }
  EXPECT_EQ(ctor_count, 1);
  EXPECT_EQ(dtor_count, 1);
}

TEST(BoxTest, UnsafeMakeUninitialized) {
  int ctor_count = 0;
  int dtor_count = 0;

  {
    Box<CtorDtorCounter> box = Box<CtorDtorCounter>::UnsafeMakeUninitialized();
    EXPECT_EQ(ctor_count, 0);
    EXPECT_EQ(dtor_count, 0);
    new (box.get()) CtorDtorCounter(ctor_count, dtor_count);

    EXPECT_EQ(ctor_count, 1);
    EXPECT_EQ(dtor_count, 0);
  }
  EXPECT_EQ(ctor_count, 1);
  EXPECT_EQ(dtor_count, 1);
};

TEST(BoxTest, FromNullableChecked) {
  NullableBox<int> a = Box<int>(12345);
  Box<int> b = Box<int>::FromNullableChecked(std::move(a));
  EXPECT_THAT(b, Pointee(12345));
}

TEST(BoxDeathTest, FromNullableCheckFail) {
  NullableBox<int> null_box = nullptr;
  EXPECT_DEATH(Box<int>::FromNullableChecked(std::move(null_box)),
               "!= nullptr");
}

TEST(BoxTest, MoveConstructor) {
  Box<int> a(12345);
  EXPECT_THAT(a, Pointee(12345));

  Box<int> b(std::move(a));
  EXPECT_THAT(b, Pointee(12345));
  EXPECT_THAT(a, IsNull());
}

class FooClass {};
class FooSubclass : public FooClass {};

TEST(BoxTest, PointerCastConstructor) {
  Box<FooSubclass> foo_sub;
  FooClass* ptr = foo_sub.get();
  Box<FooClass> foo(std::move(foo_sub));
  EXPECT_EQ(foo.get(), ptr);
}

TEST(BoxTest, MoveAssignment) {
  // Must be equivalent to swap.
  Box<int> a(12345);
  Box<int> b(0xabcde);

  a = std::move(b);

  EXPECT_THAT(a, Pointee(0xabcde));
  EXPECT_THAT(b, Pointee(12345));
}

TEST(BoxTest, IntoNullable) {
  Box<int> a(12345);
  NullableBox<int> b = std::move(a).IntoNullable();
  EXPECT_THAT(a, IsNull());
  EXPECT_THAT(b, Pointee(12345));
}

TEST(BoxTest, Swap) {
  Box<int> a(12345);
  Box<int> b(0xabcde);

  a.swap(b);

  EXPECT_THAT(a, Pointee(0xabcde));
  EXPECT_THAT(b, Pointee(12345));
}

TEST(BoxTest, StdSwap) {
  Box<int> a(12345);
  Box<int> b(0xabcde);

  using std::swap;
  swap(a, b);

  EXPECT_THAT(a, Pointee(0xabcde));
  EXPECT_THAT(b, Pointee(12345));
}

TEST(BoxTest, Dereference) {
  Box<std::string> a("foo");
  EXPECT_EQ(*a, "foo");
  EXPECT_EQ(a->size(), 3);
}

TEST(BoxTest, Equality) {
  Box<int> a(12345);
  Box<int> b(0xabcde);

  EXPECT_TRUE(a == a);
  EXPECT_TRUE(b == b);
  EXPECT_FALSE(a != a);
  EXPECT_FALSE(b != b);

  EXPECT_FALSE(a == b);
  EXPECT_FALSE(b == a);
  EXPECT_TRUE(a != b);
  EXPECT_TRUE(b != a);

  EXPECT_TRUE(a == a.get());
  EXPECT_TRUE(a != b.get());
  EXPECT_FALSE(a != a.get());
  EXPECT_FALSE(a == b.get());
}

TEST(BoxTest, NullptrEquality) {
  Box<int> a(12345);
  EXPECT_TRUE(a != nullptr);
  EXPECT_FALSE(a == nullptr);

  // Move out of a. Now a is null.
  Box<int>(std::move(a));

  EXPECT_TRUE(a == nullptr);
  EXPECT_FALSE(a != nullptr);
}

TEST(BoxTest, Ordering) {
  Box<int> a(12345);
  Box<int> b(0xabcde);

  EXPECT_EQ(a < b, a.get() < b.get());
  EXPECT_EQ(a > b, a.get() > b.get());
  EXPECT_EQ(a <= b, a.get() < b.get());
  EXPECT_EQ(a >= b, a.get() > b.get());
  EXPECT_TRUE(a <= a);
  EXPECT_TRUE(a >= a);
}

// Example of constructing a Box containing a type with a private constructor.

struct PrivateCtor {
 private:
  explicit PrivateCtor(int i) : i(i) {}

 public:
  static Box<PrivateCtor> FactoryFn(int i) {
    auto box = Box<PrivateCtor>::UnsafeMakeUninitialized();
    new (box.get()) PrivateCtor(i);  // Placement new.
    return box;
  }

  int i;
};

TEST(BoxTest, PrivateConstructor) {
  Box<PrivateCtor> box = PrivateCtor::FactoryFn(2);
  EXPECT_EQ(box->i, 2);
}

TEST(NullableBoxTest, FromBox) {
  NullableBox<int> box = Box<int>(2);
  EXPECT_THAT(box, Pointee(2));
}

TEST(NullableBoxTest, FromNullptr) {
  NullableBox<int> box = nullptr;
  EXPECT_THAT(box, IsNull());
};

TEST(NullableBoxTest, MoveConstructor) {
  NullableBox<int> a = Box<int>(12345);
  NullableBox<int> b(std::move(a));
  EXPECT_THAT(b, Pointee(12345));
  EXPECT_THAT(a, IsNull());
}

TEST(NullableBoxTest, PointerCastConstructor) {
  NullableBox<FooSubclass> foo_sub = Box<FooSubclass>();
  FooClass* ptr = foo_sub.get();
  NullableBox<FooClass> foo(std::move(foo_sub));
  EXPECT_EQ(foo.get(), ptr);
}

TEST(NullableBoxTest, MoveAssignment) {
  // Must be equivalent to swap.
  NullableBox<int> a = Box<int>(12345);
  NullableBox<int> b = Box<int>(0xabcde);

  a = std::move(b);

  EXPECT_THAT(a, Pointee(0xabcde));
  EXPECT_THAT(b, Pointee(12345));
}

TEST(NullableBoxTest, IntoBoxChecked) {
  NullableBox<int> a = Box<int>(12345);
  Box<int> b = std::move(a).AssertNotNull();
  EXPECT_THAT(a, IsNull());
  EXPECT_THAT(b, Pointee(12345));
}

TEST(NullableBoxDeathTest, IntoBoxCheckFail) {
  NullableBox<int> null_box = nullptr;
  EXPECT_DEATH(std::move(null_box).AssertNotNull(), "!= nullptr");
}

TEST(NullableBoxTest, Swap) {
  NullableBox<int> a = Box<int>(12345);
  NullableBox<int> b = Box<int>(0xabcde);

  a.swap(b);

  EXPECT_THAT(a, Pointee(0xabcde));
  EXPECT_THAT(b, Pointee(12345));
}

TEST(NullableBoxTest, StdSwap) {
  NullableBox<int> a = Box<int>(12345);
  NullableBox<int> b = Box<int>(0xabcde);

  using std::swap;
  swap(a, b);

  EXPECT_THAT(a, Pointee(0xabcde));
  EXPECT_THAT(b, Pointee(12345));
}

TEST(NullableBoxTest, Dereference) {
  NullableBox<std::string> a = Box<std::string>("foo");
  EXPECT_EQ(*a, "foo");
  EXPECT_EQ(a->size(), 3);
}

TEST(NullableBoxTest, Equality) {
  NullableBox<int> a = Box<int>(12345);
  NullableBox<int> b = Box<int>(0xabcde);

  EXPECT_TRUE(a == a);
  EXPECT_TRUE(b == b);
  EXPECT_FALSE(a != a);
  EXPECT_FALSE(b != b);

  EXPECT_FALSE(a == b);
  EXPECT_FALSE(b == a);
  EXPECT_TRUE(a != b);
  EXPECT_TRUE(b != a);

  EXPECT_TRUE(a == a.get());
  EXPECT_TRUE(a != b.get());
  EXPECT_FALSE(a != a.get());
  EXPECT_FALSE(a == b.get());
}

TEST(NullableBoxTest, NullptrEquality) {
  NullableBox<int> a = Box<int>(12345);
  EXPECT_TRUE(a != nullptr);
  EXPECT_FALSE(a == nullptr);

  // Move out of a. Now a is null.
  NullableBox<int> b = nullptr;

  EXPECT_TRUE(b == nullptr);
  EXPECT_FALSE(b != nullptr);
}

TEST(NullableBoxTest, Ordering) {
  NullableBox<int> a = Box<int>(12345);
  NullableBox<int> b = Box<int>(0xabcde);

  EXPECT_EQ(a < b, a.get() < b.get());
  EXPECT_EQ(a > b, a.get() > b.get());
  EXPECT_EQ(a <= b, a.get() < b.get());
  EXPECT_EQ(a >= b, a.get() > b.get());
  EXPECT_TRUE(a <= a);
  EXPECT_TRUE(a >= a);
}

}  // namespace base
}  // namespace principia
