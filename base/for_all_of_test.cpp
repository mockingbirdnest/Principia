#include "base/for_all_of.hpp"

#include <array>
#include <iostream>
#include <tuple>
#include <utility>

#include "gtest/gtest.h"

namespace principia {
namespace base {

using namespace principia::base::_for_all_of;

class ForAllOfTest : public ::testing::Test {
 protected:
  using AnArray = std::array<float, 3>;
  using APair = std::pair<int*, char>;
  using ATuple = std::tuple<char, double, int>;
};

TEST_F(ForAllOfTest, AnArray) {
  constexpr AnArray halved = []() {
    AnArray array{42.0, 43.0, -41.0};
    for_all_of(array).loop([](auto& value) { value /= 2; });
    return array;
  }();
  static_assert(std::get<0>(halved) == 21.0);
  static_assert(std::get<1>(halved) == 21.5);
  static_assert(std::get<2>(halved) == -20.5);
}

TEST_F(ForAllOfTest, APair) {
  static std::array<int, 2> a;
  constexpr APair incremented = []() {
    APair pair{a.data(), 'y'};
    for_all_of(pair).loop([](auto& value) { ++value; });
    return pair;
  }();
  static_assert(std::get<0>(incremented) == &a[1]);
  static_assert(std::get<1>(incremented) == 'z');
}

TEST_F(ForAllOfTest, ATuple) {
  constexpr ATuple incremented = []() {
    ATuple tuple{'a', 42.0, 666};
    for_all_of(tuple).loop([](auto& value) { ++value; });
    return tuple;
  }();
  static_assert(std::get<0>(incremented) == 'b');
  static_assert(std::get<1>(incremented) == 43.0);
  static_assert(std::get<2>(incremented) == 667);
}

TEST_F(ForAllOfTest, Parallel) {
  constexpr ATuple sum = []() {
    ATuple tuple{'A', 42.0, 666};
    AnArray array{42.0, 43.0, -41.0};
    for_all_of(tuple, array)
        .loop([](auto& tuple_element, auto const array_element) {
          tuple_element = array_element;
        });
    return tuple;
  }();
  static_assert(std::get<0>(sum) == '*');
  static_assert(std::get<1>(sum) == 43.0);
  static_assert(std::get<2>(sum) == -41);
}

TEST_F(ForAllOfTest, Indexed) {
  // In C++26 this could be a structured binding.
  constexpr auto tuple_and_array = []() {
    ATuple tuple{'A', 42.0, 666};
    AnArray array{42.0, 43.0, -41.0};
    for_all_of(tuple, array)
        .loop_indexed([]<int i>(auto& tuple_element, auto& array_element) {
          tuple_element = 64 + i;
          array_element = i;
        });
    return std::pair{tuple, array};
  }();
  constexpr auto tuple = std::get<0>(tuple_and_array);
  constexpr auto array = std::get<1>(tuple_and_array);
  static_assert(std::get<0>(tuple) == '@');
  static_assert(std::get<1>(tuple) == 65.0);
  static_assert(std::get<2>(tuple) == 66);
  static_assert(std::get<0>(array) == 0.0);
  static_assert(std::get<1>(array) == 1.0);
  static_assert(std::get<2>(array) == 2.0);
}

TEST_F(ForAllOfTest, IntegerRange) {
  // In C++26 this could be a structured binding.
  constexpr auto tuple_and_array = []() {
    ATuple tuple{'A', 42.0, 666};
    AnArray array{42.0, 43.0, -41.0};
    for_integer_range<0, 3>::loop([&]<int i> {
      std::get<i>(tuple) = 64 + i;
      std::get<i>(array) = i;
    });
    return std::pair{tuple, array};
  }();
  constexpr auto tuple = std::get<0>(tuple_and_array);
  constexpr auto array = std::get<1>(tuple_and_array);
  static_assert(std::get<0>(tuple) == '@');
  static_assert(std::get<1>(tuple) == 65.0);
  static_assert(std::get<2>(tuple) == 66);
  static_assert(std::get<0>(array) == 0.0);
  static_assert(std::get<1>(array) == 1.0);
  static_assert(std::get<2>(array) == 2.0);
}

TEST_F(ForAllOfTest, Example) {
  std::tuple const t{"a", 2.5, 3};
  std::array const a{4, 5, 6};
  for_all_of(t, a).loop([](auto const tuple_element,
                           int const array_element) {
    std::cout << tuple_element << " " << array_element << "\n";
  });
  std::tuple out{std::array{0, 0}, std::pair{"", 0}, std::pair{0.0, 0}, ""};
  for_all_of(t, a).loop_indexed([&]<int i>(auto const tuple_element,
                                           int const array_element) {
    get<(i + 1) % 3>(out) = {tuple_element, array_element};
  });
}

TEST_F(ForAllOfTest, IntegerRangeExample) {
  std::tuple t{std::string("a"), 2.5, 3};
  for_integer_range<0, 3>::loop([&]<int i> {
    if constexpr (i == 0) {
      get<i>(t) += std::to_string(i);
    } else {
      get<i>(t) += i;
    }
  });
}

}  // namespace base
}  // namespace principia
