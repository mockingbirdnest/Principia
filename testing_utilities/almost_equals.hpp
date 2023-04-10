#pragma once

#include <cfloat>
#include <cstdint>

#include <complex>
#include <string>

#include "geometry/complexification.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/point.hpp"
#include "geometry/quaternion.hpp"
#include "geometry/rotation.hpp"
#include "geometry/r3_element.hpp"
#include "geometry/r3x3_matrix.hpp"
#include "gmock/gmock.h"
#include "numerics/fixed_arrays.hpp"
#include "numerics/unbounded_arrays.hpp"

namespace principia {
namespace testing_utilities {
namespace _almost_equals {
namespace internal {

using namespace principia::quantities::_quantities;

template<typename T>
class AlmostEqualsMatcher;

// The 2-argument version of |AlmostEquals()| should always be preferred as it
// guarantees that the error bound is tight.
template<typename T>
testing::PolymorphicMatcher<AlmostEqualsMatcher<T>> AlmostEquals(
    T const& expected,
    std::int64_t max_ulps);

// The 3-argument version of |AlmostEquals()| is exclusively for use when a
// given assertion may have different errors, e.g., because it's in a loop.  It
// doesn't guarantee that the error bound is tight.  For vectors, it applies
// only to the component with the largest error.
template<typename T>
testing::PolymorphicMatcher<AlmostEqualsMatcher<T>> AlmostEquals(
    T const& expected,
    std::int64_t min_ulps,
    std::int64_t max_ulps);

template<typename T>
class AlmostEqualsMatcher final {
 public:
  explicit AlmostEqualsMatcher(T const& expected,
                               std::int64_t min_ulps,
                               std::int64_t max_ulps);

  template<typename Dimensions>
  bool MatchAndExplain(Quantity<Dimensions> const& actual,
                       testing::MatchResultListener* listener) const;
  bool MatchAndExplain(double actual,
                       testing::MatchResultListener* listener) const;
  bool MatchAndExplain(geometry::Complexification<double> actual,
                       testing::MatchResultListener* listener) const;
  template<typename Scalar>
  bool MatchAndExplain(geometry::R3Element<Scalar> const& actual,
                       testing::MatchResultListener* listener) const;
  template<typename Scalar>
  bool MatchAndExplain(geometry::R3x3Matrix<Scalar> const& actual,
                       testing::MatchResultListener* listener) const;
  bool MatchAndExplain(geometry::Quaternion const& actual,
                       testing::MatchResultListener* listener) const;
  template<typename FromFrame, typename ToFrame>
  bool MatchAndExplain(geometry::Rotation<FromFrame, ToFrame> const& actual,
                       testing::MatchResultListener* listener) const;
  template<typename Scalar, typename Frame>
  bool MatchAndExplain(geometry::Vector<Scalar, Frame> const& actual,
                       testing::MatchResultListener* listener) const;
  template<typename Scalar, typename Frame>
  bool MatchAndExplain(geometry::Bivector<Scalar, Frame> const& actual,
                       testing::MatchResultListener* listener) const;
  template<typename Scalar, typename Frame>
  bool MatchAndExplain(geometry::Trivector<Scalar, Frame> const& actual,
                       testing::MatchResultListener* listener) const;
  template<typename Vector>
  bool MatchAndExplain(geometry::Point<Vector> const& actual,
                       testing::MatchResultListener* listener) const;
  template<typename Scalar, int size>
  bool MatchAndExplain(numerics::FixedVector<Scalar, size> const& actual,
                       testing::MatchResultListener* listener) const;
  template<typename Scalar, int rows, int columns>
  bool MatchAndExplain(
      numerics::FixedMatrix<Scalar, rows, columns> const& actual,
      testing::MatchResultListener* listener) const;
  template<typename Scalar, int rows>
  bool MatchAndExplain(
      numerics::FixedLowerTriangularMatrix<Scalar, rows> const& actual,
      testing::MatchResultListener* listener) const;
  template<typename Scalar, int columns>
  bool MatchAndExplain(
      numerics::FixedUpperTriangularMatrix<Scalar, columns> const& actual,
      testing::MatchResultListener* listener) const;
  template<typename Scalar>
  bool MatchAndExplain(numerics::UnboundedVector<Scalar> const& actual,
                       testing::MatchResultListener* listener) const;
  template<typename Scalar>
  bool MatchAndExplain(
      numerics::UnboundedMatrix<Scalar> const& actual,
      testing::MatchResultListener* listener) const;
  template<typename Scalar>
  bool MatchAndExplain(
      numerics::UnboundedLowerTriangularMatrix<Scalar> const& actual,
      testing::MatchResultListener* listener) const;
  template<typename Scalar>
  bool MatchAndExplain(
      numerics::UnboundedUpperTriangularMatrix<Scalar> const& actual,
      testing::MatchResultListener* listener) const;

  void DescribeTo(std::ostream* out) const;
  void DescribeNegationTo(std::ostream* out) const;

 private:
  bool MatchAndExplainIdentical(testing::MatchResultListener* listener) const;

  T const expected_;
  std::int64_t const min_ulps_;
  std::int64_t const max_ulps_;
};

}  // namespace internal

using internal::AlmostEquals;

}  // namespace _almost_equals
}  // namespace testing_utilities
}  // namespace principia

namespace principia::testing_utilities {
using namespace principia::testing_utilities::_almost_equals;
}  // namespace principia::testing_utilities

#include "testing_utilities/almost_equals_body.hpp"
