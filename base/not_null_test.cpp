
#include "base/not_null.hpp"

#include <memory>
#include <utility>
#include <string>
#include <vector>

#include "base/macros.hpp"
#include "gtest/gtest.h"
#include "gmock/gmock.h"

using testing::Eq;

namespace principia {
namespace base {

class NotNullTest : public testing::Test {
 protected:
  // A very convoluted wrapper for the x86 add...
  void Add(not_null<int*> const destination,
           not_null<int const*> const source) {
    *destination += *source;
  }
  void Sub(not_null<int*> const destination,
           not_null<int const*> const source) {
    *destination -= *source;
  }
};

using NotNullDeathTest = NotNullTest;

TEST_F(NotNullDeathTest, DeathByNullptr) {
  EXPECT_DEATH({
    int* const null_int_ptr = nullptr;
    check_not_null(null_int_ptr);
  }, "Check failed: .* != nullptr");
  EXPECT_DEATH({
    std::unique_ptr<int> null_int_ptr;
    check_not_null(std::move(null_int_ptr));
  }, "Check failed: .* != nullptr");
}

TEST_F(NotNullTest, Move) {
  not_null<std::unique_ptr<int>> int_ptr1 = make_not_null_unique<int>(3);
#if !PRINCIPIA_COMPILER_MSVC || \
    !(_MSC_FULL_VER == 190'023'506 || \
      _MSC_FULL_VER == 190'023'918 || \
      _MSC_FULL_VER == 190'024'210 || \
      _MSC_FULL_VER == 190'024'213 || \
      _MSC_FULL_VER == 190'024'215 || \
      _MSC_FULL_VER == 191'326'215 || \
      _MSC_FULL_VER == 191'426'316 || \
      _MSC_FULL_VER == 191'426'412 || \
      _MSC_FULL_VER == 191'426'428 || \
      _MSC_FULL_VER == 191'426'429 || \
      _MSC_FULL_VER == 191'526'608 || \
      _MSC_FULL_VER == 191'526'731 || \
      _MSC_FULL_VER == 191'627'024 || \
      _MSC_FULL_VER == 191'627'025 || \
      _MSC_FULL_VER == 191'627'027 || \
      _MSC_FULL_VER == 192'027'508 || \
      _MSC_FULL_VER == 192'227'706 || \
      _MSC_FULL_VER == 192'227'724 || \
      _MSC_FULL_VER == 192'227'905 || \
      _MSC_FULL_VER == 192'328'106 || \
      _MSC_FULL_VER == 192'428'314)
  EXPECT_THAT(*(std::unique_ptr<int> const&)int_ptr1, Eq(3));
#endif
  not_null<std::unique_ptr<int>> int_ptr2 = std::move(int_ptr1);
  EXPECT_THAT(*int_ptr2, Eq(3));
  not_null<std::unique_ptr<int>> int_ptr3(std::move(int_ptr2));
  EXPECT_THAT(*int_ptr3, Eq(3));
  int_ptr2 = make_not_null_unique<int>(5);
  EXPECT_THAT(*int_ptr2, Eq(5));
  int_ptr2 = std::move(int_ptr3);
  EXPECT_THAT(*int_ptr2, Eq(3));
}

TEST_F(NotNullTest, Copy) {
  std::unique_ptr<int> const owner_of_three = std::make_unique<int>(3);
  std::unique_ptr<int> const owner_of_five = std::make_unique<int>(5);
  not_null<int*> const int_ptr1 = check_not_null(owner_of_three.get());
  not_null<int*> int_ptr2 = int_ptr1;
  not_null<int*> const int_ptr3(int_ptr2);
  EXPECT_THAT(int_ptr1, Eq(owner_of_three.get()));
  EXPECT_THAT(int_ptr2, Eq(owner_of_three.get()));
  EXPECT_THAT(int_ptr3, Eq(owner_of_three.get()));
  EXPECT_THAT(*int_ptr1, Eq(3));
  EXPECT_THAT(*int_ptr2, Eq(3));
  EXPECT_THAT(*int_ptr3, Eq(3));
  int_ptr2 = check_not_null(owner_of_five.get());
  EXPECT_THAT(*int_ptr2, Eq(5));
  int_ptr2 = int_ptr3;
  EXPECT_THAT(*int_ptr2, Eq(3));
  EXPECT_THAT(*int_ptr3, Eq(3));
}

TEST_F(NotNullTest, CheckNotNull) {
#if 0
  std::unique_ptr<int> owner_int = std::make_unique<int>(3);
  int* const constant_access_int = owner_int.get();
  not_null<int*> const constant_not_null_access_int =
      check_not_null(constant_access_int);
  check_not_null(constant_not_null_access_int);
  not_null<std::unique_ptr<int>> not_null_owner_int =
      check_not_null(std::move(owner_int));
  check_not_null(std::move(not_null_owner_int));
#endif
}

TEST_F(NotNullTest, Booleans) {
  not_null<std::unique_ptr<int>> const pointer = make_not_null_unique<int>(3);
  EXPECT_TRUE(pointer);
  EXPECT_TRUE(pointer != nullptr);
  EXPECT_FALSE(pointer == nullptr);
}

TEST_F(NotNullTest, ImplicitConversions) {
  not_null<std::unique_ptr<int>> not_null_owner_int =
      make_not_null_unique<int>(3);
  not_null<int const*> const constant_not_null_access_int =
      not_null_owner_int.get();
  // Copy constructor.
  not_null<int const*> not_null_access_constant_int =
      constant_not_null_access_int;
  // Copy assignment.
  not_null_access_constant_int = not_null_owner_int.get();
  // Move constructor.
  not_null<std::unique_ptr<int const>> not_null_owner_constant_int =
      make_not_null_unique<int>(5);
  // Move assignment.
  not_null_owner_constant_int = std::move(not_null_owner_int);
}

TEST_F(NotNullTest, Arrow) {
  not_null<std::unique_ptr<std::string>> not_null_owner_string =
      make_not_null_unique<std::string>("-");
  not_null_owner_string->append(">");
  EXPECT_THAT(*not_null_owner_string, Eq("->"));
  not_null<std::string*> not_null_access_string = not_null_owner_string.get();
  not_null_access_string->insert(0, "operator");
  EXPECT_THAT(*not_null_access_string, Eq("operator->"));
}

TEST_F(NotNullTest, NotNullNotNull) {
  not_null<std::unique_ptr<int>> owner = make_not_null_unique<int>(42);
  not_null<not_null<int*>> x = owner.get();
  not_null<int*> y = x;
  *x = 6 * 9;
  EXPECT_THAT(*owner, Eq(6 * 9));
  EXPECT_THAT(*x, Eq(6 * 9));
  EXPECT_THAT(*y, Eq(6 * 9));
}

TEST_F(NotNullTest, CheckArguments) {
  not_null<std::unique_ptr<int>> accumulator = make_not_null_unique<int>(21);
  std::unique_ptr<int const> twenty_one = std::make_unique<int const>(21);
  Add(accumulator.get(), check_not_null(twenty_one.get()));
  EXPECT_THAT(*twenty_one, Eq(21));
  EXPECT_THAT(*accumulator, Eq(42));
#if 0
  Sub(check_not_null(accumulator.get()), check_not_null(twenty_one.get()));
  EXPECT_THAT(*accumulator, Eq(21));
#endif
}

TEST_F(NotNullTest, DynamicCast) {
  class Base {
   public:
    virtual ~Base() = default;
  };
  class Derived : public Base {};
  {
    not_null<Base*> b = new Derived;
    [[maybe_unused]] not_null<Derived*> d = dynamic_cast_not_null<Derived*>(b);
    delete b;
  }
  {
    not_null<std::unique_ptr<Base>> b = make_not_null_unique<Derived>();
    not_null<std::unique_ptr<Derived>> d =
        dynamic_cast_not_null<std::unique_ptr<Derived>>(std::move(b));
  }
}

TEST_F(NotNullTest, RValue) {
  std::unique_ptr<int> owner_int = std::make_unique<int>(1);
  not_null<int*> not_null_int = new int(2);
  EXPECT_NE(owner_int.get(), not_null_int);

  not_null<std::unique_ptr<int>> not_null_owner_int =
      make_not_null_unique<int>(4);
  std::vector<int*> v1;
  // |v1.push_back| would be ambiguous here if |not_null<pointer>| had an
  // |operator pointer const&() const&| instead of an
  // |operator pointer const&&() const&|.
  v1.push_back(not_null_owner_int.get());
  // |emplace_back| is fine no matter what.
  v1.emplace_back(not_null_owner_int.get());
  EXPECT_EQ(4, *v1[0]);
  EXPECT_EQ(4, *not_null_owner_int);

  std::vector<int*> v2;
  // NOTE(egg): The following fails using clang, I'm not sure where the bug is.
  // More generally, if functions |foo(int*&&)| and |foo(int* const&)| exist,
  // |foo(non_rvalue_not_null)| fails to compile with clang ("no viable
  // conversion"), but without the |foo(int*&&)| overload it compiles.
  // This is easily circumvented using a temporary (whereas the ambiguity that
  // would result from having an |operator pointer const&() const&| would
  // entirely prevent conversion of |not_null<unique_ptr<T>>| to
  // |unique_ptr<T>|), so we ignore it.
#if PRINCIPIA_COMPILER_MSVC
  v2.push_back(not_null_int);
#else
  int* const temporary = not_null_int;
  v2.push_back(temporary);
#endif
  EXPECT_EQ(2, *v2[0]);
  EXPECT_EQ(2, *not_null_int);

  // |std::unique_ptr<int>::operator=| would be ambiguous here between the move
  // and copy assignments if |not_null<pointer>| had an
  // |operator pointer const&() const&| instead of an
  // |operator pointer const&&() const&|.
  // Note that one of the overloads (the implicit copy assignment operator) is
  // deleted, but this does not matter for overload resolution; however MSVC
  // compiles this even when it is ambiguous, while clang correctly fails.
  owner_int = make_not_null_unique<int>(1729);

  not_null<int*> not_null_access_int = owner_int.get();
  [[maybe_unused]] int const* access_const_int = not_null_access_int;

  not_null<std::shared_ptr<int>> not_null_shared_int = std::make_shared<int>(3);
  // This exercises |operator OtherPointer() const&|.  The conversion from
  // |not_null<int*>| to |int const*| does not; instead it goes through the
  // conversion |not_null<int*>| -> |int*| -> |int const*|, where the latter is
  // a qualification adjustment, part of a standard conversion sequence.
  std::shared_ptr<int const> shared_const_int = not_null_shared_int;

  // MSVC seems to be confused by templatized move-conversion operators.
#if PRINCIPIA_COMPILER_MSVC
  // Uses |not_null(not_null<OtherPointer>&& other)| followed by
  // |operator pointer&&() &&|, so the first conversion has to be explicit.
  std::unique_ptr<int const> owner_const_int =
      static_cast<not_null<std::unique_ptr<int const>>>(
          make_not_null_unique<int>(5));
#else
  // Uses |operator OtherPointer() &&|.
  std::unique_ptr<int const> owner_const_int = make_not_null_unique<int>(5);
#endif
  EXPECT_EQ(5, *owner_const_int);
}

}  // namespace base
}  // namespace principia
