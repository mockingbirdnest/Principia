#pragma once

#include <concepts>
#include <utility>

#include "numerics/concepts.hpp"
#include "numerics/transposed_view.hpp"
#include "quantities/named_quantities.hpp"

namespace principia {
namespace numerics {
namespace _matrix_views {
namespace internal {

// TODO(phl): The view stuff should be (1) made complete, i.e., have all the
// operations that exist for fixed/unbounded vectors/matrices; (2) unified with
// fixed/unbounded arrays so that we don't have to write each algorithm N times;
// (3) tested.

using namespace principia::numerics::_concepts;
using namespace principia::numerics::_transposed_view;
using namespace principia::quantities::_named_quantities;

// A view of a rectangular block of a matrix.  This view is |two_dimensional|.
template<typename Matrix>
  requires two_dimensional<Matrix>
struct BlockView {
  using Scalar = typename Matrix::Scalar;

  Matrix& matrix;
  int first_row;
  int last_row;
  int first_column;
  int last_column;

  constexpr Scalar& operator()(int row, int column);
  constexpr Scalar const& operator()(int row, int column) const;

  template<typename T>
    requires two_dimensional<T> &&
             std::same_as<typename T::Scalar, typename Matrix::Scalar>
  BlockView& operator-=(T const& right);

  constexpr int rows() const;
  constexpr int columns() const;
};

// A view of a column of a matrix.  This view is |one_dimensional|.
template<typename Matrix>
  requires two_dimensional<Matrix>
struct ColumnView {
  using Scalar = typename Matrix::Scalar;

  Matrix& matrix;
  int first_row;
  int last_row;
  int column;

  constexpr Scalar& operator[](int index);
  constexpr Scalar const& operator[](int index) const;

  ColumnView& operator/=(double right);

  Scalar Norm() const;
  Square<Scalar> Norm²() const;
  constexpr int size() const;
};

template<typename Matrix>
std::ostream& operator<<(std::ostream& out,
                         ColumnView<Matrix> const& view);

}  // namespace internal

using internal::BlockView;
using internal::ColumnView;

}  // namespace _matrix_views
}  // namespace numerics
}  // namespace principia

#include "numerics/matrix_views_body.hpp"