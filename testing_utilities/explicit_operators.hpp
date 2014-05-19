// This file contains template functions used to refer unambiguously to
// overloaded operators, for instance when function pointers are needed.

#pragma once

namespace principia {
namespace testing_utilities {

template<typename ResultType, typename LeftType, typename RightType>
ResultType Plus(LeftType const& left, RightType const& right);

template<typename ResultType, typename LeftType, typename RightType>
ResultType Minus(LeftType const& left, RightType const& right);

template<typename ResultType, typename LeftType, typename RightType>
ResultType Times(LeftType const& left, RightType const& right);

template<typename ResultType, typename LeftType, typename RightType>
ResultType Divide(LeftType const& left, RightType const& right);

}  // namespace testing_utilities
}  // namespace principia

#include "TestUtilities/ExplicitOperators-body.hpp"
