#pragma once

namespace principia {
namespace testing_utilities {

template<typename ResultType, typename RightType>
inline ResultType Minus(RightType const& right) {
  return -right;
}

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

}  // namespace testing_utilities
}  // namespace principia
