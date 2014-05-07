#pragma once

namespace principia {
namespace test_utilities {

template<typename ResultType, typename LeftType, typename RightType>
inline ResultType Plus(LeftType const& left, RightType const& right) {
  return left + right;
}

template<typename ResultType, typename LeftType, typename RightType>
inline ResultType Minus(LeftType const& left, RightType const& right) {
  return left - right;
}

template<typename ResultType, typename LeftType, typename RightType>
inline ResultType Times(LeftType const& left, RightType const& right) {
  return left * right;
}

template<typename ResultType, typename LeftType, typename RightType>
inline ResultType Divide(LeftType const& left, RightType const& right) {
  return left / right;
}

}  // test_utilities
}  // principia
