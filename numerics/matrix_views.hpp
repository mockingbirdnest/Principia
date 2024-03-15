#pragma once

#include "numerics/concepts.hpp"
#include "numerics/transposed_view.hpp"
#include "numerics/unbounded_arrays.hpp"
#include "quantities/named_quantities.hpp"

namespace principia {
namespace numerics {
namespace _matrix_views {
namespace internal {

// TODO(phl): The view stuff should be (1) made completed, i.e., have all the
// operations that exist for fixed/unbounded vectors/matrices; (2) unified with
// fixed/unbounded arrays so that we don't have to write each algorithm N times.

using namespace principia::numerics::_concepts;
using namespace principia::numerics::_transposed_view;
using namespace principia::numerics::_unbounded_arrays;
using namespace principia::quantities::_named_quantities;

// A view of a column of a matrix.  This view is `two_dimensional`.
template<typename Matrix>
  requires two_dimensional<Matrix>
struct BlockView {
  using Scalar = typename Matrix::Scalar;

  Matrix& matrix;
  int first_row;
  int last_row;
  int first_column;
  int last_column;

  constexpr int rows() const;
  constexpr int columns() const;

  constexpr Scalar& operator()(int row, int column);
  constexpr Scalar const& operator()(int row, int column) const;

  BlockView& operator-=(UnboundedMatrix<Scalar> const& right);
};

// A view of a column of a matrix.  This view is `one_dimensional`.
template<typename Matrix>
  requires two_dimensional<Matrix>
struct ColumnView {
  using Scalar = typename Matrix::Scalar;

  Matrix& matrix;
  int first_row;
  int last_row;
  int column;

  Scalar Norm() const;
  Square<Scalar> Norm²() const;
  constexpr int size() const;

  // Constructs an unbounded vector by copying data from the view.  Note that
  // the result is unbounded even if the matrix being viewed is a FixedMatrix.
  //TODO(phl)Move
  explicit operator UnboundedVector<Scalar>() const;

  constexpr Scalar& operator[](int index);
  constexpr Scalar const& operator[](int index) const;

  ColumnView& operator/=(double right);
};

//TODO(phl)Move these operators.

template<typename LMatrix, typename RScalar>
  requires two_dimensional<LMatrix>
UnboundedVector<Product<typename LMatrix::Scalar, RScalar>> operator*(
    BlockView<LMatrix> const& left,
    UnboundedVector<RScalar> const& right);

template<typename LMatrix, typename RScalar>
  requires two_dimensional<LMatrix>
UnboundedVector<Product<typename LMatrix::Scalar, RScalar>> operator*(
    TransposedView<BlockView<LMatrix>> const& left,
    UnboundedVector<RScalar> const& right);

}  // namespace internal

using internal::BlockView;
using internal::ColumnView;

}  // namespace _matrix_views
}  // namespace numerics
}  // namespace principia

#include "numerics/matrix_views_body.hpp"
