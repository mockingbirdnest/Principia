#pragma once

namespace principia {
namespace test_utilities {

template<typename ResultType, typename LeftType, typename RightType>
ResultType Plus(LeftType const& left, RightType const& right);

template<typename ResultType, typename LeftType, typename RightType>
ResultType Minus(LeftType const& left, RightType const& right);

template<typename ResultType, typename LeftType, typename RightType>
ResultType Times(LeftType const& left, RightType const& right);

template<typename ResultType, typename LeftType, typename RightType>
ResultType Divide(LeftType const& left, RightType const& right);

}  // test_utilities
}  // principia

#include "TestUtilities/ExplicitOperators-body.hpp"
