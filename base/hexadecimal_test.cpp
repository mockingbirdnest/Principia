#pragma once

#include "base/hexadecimal.hpp"

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "glog/stl_logging.h"

namespace principia {
namespace base {

class HexadecimalTest : public testing::Test {
 protected:
  std::string const bytes_ = std::string("\0\x7F\x80\xFFgh\n\7", 8);
  std::string const lowercase_digits_ = "00""7f""80""ff""67""68""0a""07";
  std::string const uppercase_digits_ = "00""7F""80""FF""67""68""0A""07";
};

TEST_F(HexadecimalTest, EncodeAndDecode) {
  std::string digits;
  HexadecimalEncode<std::string>(bytes_, &digits);
  EXPECT_EQ(uppercase_digits_, digits);
  std::string bytes;
  HexadecimalDecode<std::string>(digits, &bytes);
  EXPECT_EQ(bytes_, bytes);
}

TEST_F(HexadecimalTest, InPlace) {
  std::string str = bytes_;
  HexadecimalEncode<std::string>(str, &str);
  EXPECT_EQ(uppercase_digits_, str);
  HexadecimalDecode<std::string>(str, &str);
  EXPECT_EQ(bytes_, str);
}

TEST_F(HexadecimalTest, CaseInsensitive) {
  std::string bytes;
  HexadecimalDecode<std::string>(lowercase_digits_, &bytes);
  EXPECT_EQ(bytes_, bytes);
  HexadecimalDecode<std::string>(uppercase_digits_, &bytes);
  EXPECT_EQ(bytes_, bytes);
}

TEST_F(HexadecimalTest, Validity) {
  std::string bytes;
  std::string digits = uppercase_digits_;
  bool valid;
  valid = HexadecimalDecode<std::string>(digits, &bytes);
  EXPECT_TRUE(valid);
  EXPECT_EQ(bytes_, bytes);
  digits = "abc";
  valid = HexadecimalDecode<std::string>(digits, &bytes);
  EXPECT_FALSE(valid);
  EXPECT_EQ("\xAB", bytes);
  digits = "0agcde";
  valid = HexadecimalDecode<std::string>(digits, &bytes);
  EXPECT_FALSE(valid);
  EXPECT_EQ(std::string("\x0A\x00\xDE", 3), bytes);
}

// Allocator adaptor that interposes construct() calls to
// convert value initialization into default initialization.
template <typename T, typename A=std::allocator<T>>
class default_init_allocator : public A {
  typedef std::allocator_traits<A> a_t;
public:
  template <typename U> struct rebind {
    using other =
      default_init_allocator<
        U, typename a_t::template rebind_alloc<U>
      >;
  };

  template <typename U>
  void construct(U* ptr) {
    ::new (static_cast<void*>(ptr)) U;
  }
  template <typename U, typename...Args>
  void construct(U* ptr, Args&&... args) {
    a_t::construct(static_cast<A&>(*this),
                   ptr, std::forward<Args>(args)...);
  }
};

template<typename T>
using uninitialized_vector = std::vector<T, default_init_allocator<T>>;

TEST_F(HexadecimalTest, Foo) {
    {
  uninitialized_vector<uint8_t> foo;
  foo.reserve(10000);
  memset(foo.data(), 42, 1000);
  foo.resize(10000);
  std::vector<uint8_t> bar;
  bar.reserve(1000);
  memset(bar.data(), 42, 1000);
  bar.resize(1000);
  LOG(ERROR)<<foo;
  LOG(ERROR)<<bar;
    }

  std::clock_t t_uninitialized = 0;
  uint8_t n = 42;
  std::size_t const size = 1 << 25;
  std::size_t const iterations = 100;
  unsigned long total = 0;
  for (int i = iterations; i > 0; --i) {
    t_uninitialized -= std::clock();
    uninitialized_vector<uint8_t> foo(size);
    for (auto& x : foo) {
      x = ++n;
    }
    t_uninitialized += std::clock();
    for (auto const x : foo) {
      total += x;
    }
  }
  LOG(ERROR)<<t_uninitialized;
  std::clock_t t_initialized = 0;
  for (int i = iterations; i > 0; --i) {
    t_initialized -= std::clock();
    std::vector<uint8_t> foo(size);
    for (auto& x : foo) {
      x = ++n;
    }
    t_initialized += std::clock();
    for (auto const x : foo) {
      total += x;
    }
  }
  LOG(ERROR)<<t_initialized;
  std::clock_t t_reserved = 0;
  for (int i = iterations; i > 0; --i) {
    t_reserved -= std::clock();
    std::vector<uint8_t> foo;
    foo.reserve(size);
    for (int j = 0; j < size; ++j) {
      foo.push_back(++n);
    }
    t_reserved += std::clock();
    for (auto const x : foo) {
      total += x;
    }
  }
  LOG(ERROR)<<t_reserved;
  LOG(ERROR)<<total;
}

TEST_F(HexadecimalTest, Bar) {;
  std::size_t const size = 1 << 25;
  std::size_t const iterations = 100;
  char const* const hi = "hi";
  std::vector<uint8_t> greet_a_lot(size);
  for (uint8_t* p = &greet_a_lot[0]; p < &greet_a_lot[size]; ++p) {
      *p = *hi;
      *++p = *(hi + 1);
    };
  LOG(ERROR)<<greet_a_lot.front()<<"..."<<greet_a_lot.back();
  std::clock_t t_memcpy = 0;
  for (int i = iterations; i > 0; --i) {
    std::vector<uint8_t> foo(size);
    t_memcpy -= std::clock();
    for (uint8_t* p = &foo[0]; p < &foo[size]; p += 2) {
      memcpy(p, hi, 2);
    }
    t_memcpy += std::clock();
    EXPECT_EQ(greet_a_lot, foo);
  }
  LOG(ERROR)<<t_memcpy;
  std::clock_t t_char = 0;
  for (int i = iterations; i > 0; --i) {
    std::vector<uint8_t> foo(size);
    t_char -= std::clock();
    for (uint8_t* p = &foo[0]; p < &foo[size]; ++p) {
      *p = *hi;
      *++p = *(hi + 1);
    }
    t_char += std::clock();
    EXPECT_EQ(greet_a_lot, foo);
  }
  LOG(ERROR)<<t_char;
}

}  // namespace base
}  // namespace principia
