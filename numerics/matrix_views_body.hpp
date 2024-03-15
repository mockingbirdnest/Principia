#pragma once

#include "numerics/matrix_views.hpp"

#include "base/tags.hpp"
#include "quantities/elementary_functions.hpp"

namespace principia {
namespace numerics {
namespace _matrix_views {
namespace internal {

using namespace principia::base::_tags;
using namespace principia::quantities::_elementary_functions;

template<typename Matrix>
  requires two_dimensional<Matrix>
constexpr auto BlockView<Matrix>::rows() const -> int {
  return last_row - first_row + 1;
}

template<typename Matrix>
  requires two_dimensional<Matrix>
constexpr auto BlockView<Matrix>::columns() const -> int{
  return last_column - first_column + 1;
}

template<typename Matrix>
  requires two_dimensional<Matrix>
constexpr auto BlockView<Matrix>::operator()(
    int const row, int const column) -> Scalar& {
  CONSTEXPR_DCHECK(row <= last_row - first_row);
  CONSTEXPR_DCHECK(column <= last_column - first_column);
  return matrix(first_row + row, first_column + column);
}

template<typename Matrix>
  requires two_dimensional<Matrix>
constexpr auto BlockView<Matrix>::operator()(
    int const row,
    int const column) const -> Scalar const& {
  CONSTEXPR_DCHECK(row <= last_row - first_row);
  CONSTEXPR_DCHECK(column <= last_column - first_column);
  return matrix(first_row + row, first_column + column);
}

template<typename Matrix>
  requires two_dimensional<Matrix>
auto BlockView<Matrix>::operator-=(
    UnboundedMatrix<Scalar> const& right) -> BlockView<Matrix>& {
  CHECK_EQ(rows(), right.rows());
  CHECK_EQ(columns(), right.columns());
  for (int i = 0; i < right.rows(); ++i) {
    for (int j = 0; j < right.columns(); ++j) {
      matrix(first_row + i, first_column + j) -= right(i, j);
    }
  }
  return *this;
}


template<typename Matrix>
  requires two_dimensional<Matrix>
auto ColumnView<Matrix>::Norm() const -> Scalar {
  return Sqrt(Norm²());
}

template<typename Matrix>
  requires two_dimensional<Matrix>
auto ColumnView<Matrix>::Norm²() const -> Square<Scalar> {
  Square<Scalar> result{};
  for (int i = first_row; i <= last_row; ++i) {
    result += Pow<2>(matrix(i, column));
  }
  return result;
}

template<typename Matrix>
  requires two_dimensional<Matrix>
constexpr auto ColumnView<Matrix>::size() const -> int {
  return last_row - first_row + 1;
}

template<typename Matrix>
  requires two_dimensional<Matrix>
ColumnView<Matrix>::operator UnboundedVector<typename Matrix::Scalar>() const {
  UnboundedVector<Scalar> result(size(), uninitialized);
  for (int i = first_row; i <= last_row; ++i) {
    result[i - first_row] = matrix(i, column);
  }
  return result;
}

template<typename Matrix>
  requires two_dimensional<Matrix>
constexpr auto ColumnView<Matrix>::operator[](int const index) -> Scalar& {
  CONSTEXPR_DCHECK(index <= last_row - first_row);
  return matrix(first_row + index, column);
}

template<typename Matrix>
  requires two_dimensional<Matrix>
constexpr auto ColumnView<Matrix>::operator[](
    int const index) const -> Scalar const& {
  CONSTEXPR_DCHECK(index <= last_row - first_row);
  return matrix(first_row + index, column);
}

template<typename Matrix>
  requires two_dimensional<Matrix>
auto ColumnView<Matrix>::operator/=(double const right) -> ColumnView<Matrix>& {
  for (int i = first_row; i < last_row; ++i) {
    matrix(i, column) /= right;
  }
}

template<typename Matrix>
  requires two_dimensional<Matrix>
UnboundedVector<double> Normalize(ColumnView<Matrix> const& view) {
  return UnboundedVector<typename Matrix::Scalar>(view) / view.Norm();
}

template<typename Matrix>
  requires two_dimensional<Matrix>
std::ostream& operator<<(std::ostream& out,
                         ColumnView<Matrix> const& view) {
  std::stringstream s;
  for (int i = 0; i < view.size(); ++i) {
    s << (i == 0 ? "{" : "") << view[i]
      << (i == view.size() - 1 ? "}" : ", ");
  }
  out << s.str();
  return out;
}


template<typename LMatrix, typename RScalar>
  requires two_dimensional<LMatrix>
UnboundedVector<Product<typename LMatrix::Scalar, RScalar>> operator*(
    BlockView<LMatrix> const& left,
    UnboundedVector<RScalar> const& right) {
  CHECK_EQ(left.columns(), right.size());
  UnboundedVector<Product<typename LMatrix::Scalar, RScalar>> result(
      left.rows());
  for (int i = 0; i < left.rows(); ++i) {
    auto& result_i = result[i];
    for (int j = 0; j < left.columns(); ++j) {
      result_i += left(i, j) * right[j];
    }
  }
  return result;
}

template<typename LMatrix, typename RScalar>
  requires two_dimensional<LMatrix>
UnboundedVector<Product<typename LMatrix::Scalar, RScalar>> operator*(
    TransposedView<BlockView<LMatrix>> const& left,
    UnboundedVector<RScalar> const& right) {
  CHECK_EQ(left.columns(), right.size());
  UnboundedVector<Product<typename LMatrix::Scalar, RScalar>> result(
      left.rows());
  for (int i = 0; i < left.rows(); ++i) {
    auto& result_i = result[i];
    for (int j = 0; j < left.columns(); ++j) {
      result_i += left(i, j) * right[j];
    }
  }
  return result;
}

}  // namespace internal
}  // namespace _matrix_views
}  // namespace numerics
}  // namespace principia
