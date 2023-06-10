#include "base/non_static_thread_local.hpp"

#include <thread>
#include <vector>

#include "gtest/gtest.h"
#include "gmock/gmock.h"

namespace principia {
namespace base {

using ::testing::Eq;
using ::testing::IsEmpty;
using ::testing::IsNull;
using ::testing::Pointee;
using ::testing::SizeIs;
using namespace principia::base::_non_static_thread_local;
using namespace std::chrono_literals;

class NonStaticThreadLocalTest : public ::testing::Test {
 protected:
  struct T {
    non_static_thread_local<int> i = 9;
    non_static_thread_local<int> j = 8;
    non_static_thread_local<std::unique_ptr<int>> p;
    non_static_thread_local<std::vector<int>> v =
        std::initializer_list{1, 2, 5};
  };

  T x_;
  T y_;
};

TEST_F(NonStaticThreadLocalTest, ReadWrite) {
  std::thread t1([this] {
    EXPECT_THAT(x_.i(), Eq(9));
    EXPECT_THAT(y_.i(), Eq(9));
    EXPECT_THAT(x_.j(), Eq(8));
    EXPECT_THAT(x_.v(), SizeIs(3));
    x_.v().push_back(3);
    x_.i() = 1;
    y_.i() = 11;
    std::this_thread::sleep_for(10ms);
    EXPECT_THAT(x_.i(), Eq(1));
    EXPECT_THAT(y_.i(), Eq(11));
    EXPECT_THAT(x_.v(), SizeIs(4));
    EXPECT_THAT(x_.j(), Eq(8));
  });
  std::thread t2([this] {
    EXPECT_THAT(x_.i(), Eq(9));
    EXPECT_THAT(x_.v(), SizeIs(3));
    x_.v().clear();
    x_.i() = 2;
    x_.j() = 102;
    y_.i() = 12;
    EXPECT_THAT(x_.i(), Eq(2));
    EXPECT_THAT(y_.i(), Eq(12));
    EXPECT_THAT(x_.v(), IsEmpty());
    EXPECT_THAT(x_.j(), Eq(102));
  });
  t1.join();
  t2.join();
}

TEST_F(NonStaticThreadLocalTest, Move) {
  std::thread t1([this] {
    x_.p() = std::make_unique<int>(1);
    y_.p() = std::move(x_).p();
    EXPECT_THAT(y_.p(), Pointee(1));
    EXPECT_THAT(x_.p(), IsNull());
  });
  std::thread t2([this] {
    y_.p() = std::make_unique<int>(2);
    x_.p() = std::move(y_).p();
    EXPECT_THAT(x_.p(), Pointee(2));
    EXPECT_THAT(y_.p(), IsNull());
  });
  t1.join();
  t2.join();
}

}  // namespace base
}  // namespace principia
