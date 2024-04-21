#pragma once

#include "numerics/matrix_views.hpp"

#include <utility>

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
template<typename T>
  requires two_dimensional<T> && same_elements_as<T, Matrix>
BlockView<Matrix>& BlockView<Matrix>::operator=(T const& right) {
  DCHECK_EQ(rows(), right.rows());
  DCHECK_EQ(columns(), right.columns());
  for (int i = 0; i < rows(); ++i) {
    for (int j = 0; j < columns(); ++j) {
      matrix(first_row + i, first_column + j) = right(i, j);
    }
  }
  return *this;
}

template<typename Matrix>
  requires two_dimensional<Matrix>
template<typename T>
  requires two_dimensional<T> && same_elements_as<T, Matrix>
BlockView<Matrix>& BlockView<Matrix>::operator=(T&& right) {
  DCHECK_EQ(rows(), right.rows());
  DCHECK_EQ(columns(), right.columns());
  for (int i = 0; i < rows(); ++i) {
    for (int j = 0; j < columns(); ++j) {
      matrix(first_row + i, first_column + j) = std::move(right(i, j));
    }
  }
  return *this;
}

template<typename Matrix>
  requires two_dimensional<Matrix>
template<typename T>
  requires two_dimensional<T> && same_elements_as<T, Matrix>
BlockView<Matrix>& BlockView<Matrix>::operator+=(T const& right) {
  DCHECK_EQ(rows(), right.rows());
  DCHECK_EQ(columns(), right.columns());
  for (int i = 0; i < rows(); ++i) {
    for (int j = 0; j < columns(); ++j) {
      matrix(first_row + i, first_column + j) += right(i, j);
    }
  }
  return *this;
}

template<typename Matrix>
  requires two_dimensional<Matrix>
template<typename T>
  requires two_dimensional<T> && same_elements_as<T, Matrix>
auto BlockView<Matrix>::operator-=(T const& right) -> BlockView<Matrix>& {
  DCHECK_EQ(rows(), right.rows());
  DCHECK_EQ(columns(), right.columns());
  for (int i = 0; i < rows(); ++i) {
    for (int j = 0; j < columns(); ++j) {
      matrix(first_row + i, first_column + j) -= right(i, j);
    }
  }
  return *this;
}

template<typename Matrix>
  requires two_dimensional<Matrix>
BlockView<Matrix>& BlockView<Matrix>::operator*=(double const right) {
  for (int i = 0; i < rows(); ++i) {
    for (int j = 0; j < columns(); ++j) {
      matrix(first_row + i, first_column + j) *= right;
    }
  }
  return *this;
}

template<typename Matrix>
  requires two_dimensional<Matrix>
BlockView<Matrix>& BlockView<Matrix>::operator/=(double const right) {
  for (int i = 0; i < rows(); ++i) {
    for (int j = 0; j < columns(); ++j) {
      matrix(first_row + i, first_column + j) /= right;
    }
  }
  return *this;
}

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
template<typename T>
  requires one_dimensional<T> && same_elements_as<T, Matrix>
ColumnView<Matrix>& ColumnView<Matrix>::operator=(T const& right) {
  DCHECK_EQ(size(), right.size());
  for (int i = 0; i < size(); ++i) {
    matrix(first_row + i, column) = right[i];
  }
  return *this;
}

template<typename Matrix>
  requires two_dimensional<Matrix>
template<typename T>
  requires one_dimensional<T> && same_elements_as<T, Matrix>
ColumnView<Matrix>& ColumnView<Matrix>::operator=(T&& right) {
  DCHECK_EQ(size(), right.size());
  for (int i = 0; i < size(); ++i) {
    matrix(first_row + i, column) = std::move(right[i]);
  }
  return *this;
}

template<typename Matrix>
  requires two_dimensional<Matrix>
template<typename T>
  requires one_dimensional<T> && same_elements_as<T, Matrix>
ColumnView<Matrix>& ColumnView<Matrix>::operator+=(T const& right) {
  DCHECK_EQ(size(), right.size());
  for (int i = 0; i < size(); ++i) {
    matrix(first_row + i, column) += right[i];
  }
  return *this;
}

template<typename Matrix>
  requires two_dimensional<Matrix>
template<typename T>
  requires one_dimensional<T> && same_elements_as<T, Matrix>
ColumnView<Matrix>& ColumnView<Matrix>::operator-=(T const& right) {
  DCHECK_EQ(size(), right.size());
  for (int i = 0; i < size(); ++i) {
    matrix(first_row + i, column) -= right[i];
  }
  return *this;
}

template<typename Matrix>
  requires two_dimensional<Matrix>
ColumnView<Matrix>& ColumnView<Matrix>::operator*=(double const right) {
  for (int i = first_row; i <= last_row; ++i) {
    matrix(i, column) *= right;
  }
  return *this;
}

template<typename Matrix>
  requires two_dimensional<Matrix>
ColumnView<Matrix>& ColumnView<Matrix>::operator/=(double const right) {
  for (int i = first_row; i <= last_row; ++i) {
    matrix(i, column) /= right;
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

template<typename LMatrix, typename RMatrix>
Product<typename LMatrix::Scalar, typename RMatrix::Scalar> operator*(
    TransposedView<ColumnView<LMatrix>> const& left,
    ColumnView<RMatrix> const& right) {
  DCHECK_EQ(left.size(), right.size());
  Product<typename LMatrix::Scalar, typename RMatrix::Scalar> result{};
  for (int i = 0; i < left.size(); ++i) {
    result += left[i] * right[i];
  }
  return result;
}

template<typename Matrix, typename T>
  requires two_dimensional<T> && same_elements_as<T, Matrix>
bool operator==(BlockView<Matrix> const& left, T const& right) {
  if (left.rows() != right.rows() || left.columns() != right.columns()) {
    return false;
  }
  for (int i = 0; i < left.rows(); ++i) {
    for (int j = 0; j < left.columns(); ++j) {
      if (left(i, j) != right(i, j)) {
        return false;
      }
    }
  }
  return true;
}

template<typename Matrix, typename T>
  requires one_dimensional<T> && same_elements_as<T, Matrix>
bool operator==(ColumnView<Matrix> const& left, T const& right) {
  if (left.size() != right.size()) {
    return false;
  }
  for (int i = 0; i < left.size(); ++i) {
    if (left[i] != right[i]) {
      return false;
    }
  }
  return true;
}

template<typename Matrix>
std::ostream& operator<<(std::ostream& out, BlockView<Matrix> const& view) {
  out << "rows: " << view.rows() << " columns: " << view.columns() << "\n";
  for (int i = 0; i < view.rows(); ++i) {
    out << "{";
    for (int j = 0; j < view.columns(); ++j) {
      out << view(i, j);
      if (j < view.columns() - 1) {
        out << ", ";
      }
    }
    out << "}\n";
  }
  return out;
}

template<typename Matrix>
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

}  // namespace internal
}  // namespace _matrix_views
}  // namespace numerics
}  // namespace principia
