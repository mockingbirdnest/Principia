
#pragma once

#include "numerics/unbounded_arrays.hpp"

#include <cmath>
#include <vector>

#include "base/macros.hpp"
#include "quantities/elementary_functions.hpp"
#include "unbounded_arrays.hpp"

namespace principia {
namespace numerics {
namespace internal_unbounded_arrays {

using base::uninitialized;
using quantities::Pow;
using quantities::Sqrt;

template<class T>
template<class U, class... Args>
inline void uninitialized_allocator<T>::construct(U* const p, Args&&... args) {
  ::new ((void*)p) U(std::forward<Args>(args)...);  // NOLINT
}

template<typename Scalar>
UnboundedVector<Scalar>::UnboundedVector(int const size)
    : data_(size, Scalar{}) {}

template<typename Scalar>
UnboundedVector<Scalar>::UnboundedVector(int const size, uninitialized_t)
    : data_(size)  {}

template<typename Scalar>
UnboundedVector<Scalar>::UnboundedVector(std::initializer_list<Scalar> data)
    : data_(std::move(data)) {}

template<typename Scalar>
void UnboundedVector<Scalar>::Extend(int const extra_size) {
  DCHECK_LE(0, extra_size);
  data_.resize(data_.size() + extra_size, Scalar{});
}

template<typename Scalar>
void UnboundedVector<Scalar>::Extend(int const extra_size, uninitialized_t) {
  DCHECK_LE(0, extra_size);
  data_.resize(data_.size() + extra_size);
}

template<typename Scalar>
void UnboundedVector<Scalar>::Extend(std::initializer_list<Scalar> data) {
  std::move(data.begin(), data.end(), std::back_inserter(data_));
}

template<typename Scalar>
void UnboundedVector<Scalar>::EraseToEnd(int const begin_index) {
  data_.erase(data_.begin() + begin_index, data_.end());
}

template<typename Scalar>
int UnboundedVector<Scalar>::size() const {
  return data_.size();
}

template<typename Scalar>
bool UnboundedVector<Scalar>::operator==(UnboundedVector const& right) const {
  return data_ == right.data_;
}

template<typename Scalar>
Scalar& UnboundedVector<Scalar>::operator[](int const index) {
  return data_[index];
}

template<typename Scalar>
Scalar const& UnboundedVector<Scalar>::operator[](int const index) const {
  return data_[index];
}

template<typename Scalar>
UnboundedLowerTriangularMatrix<Scalar>::UnboundedLowerTriangularMatrix(
    int const rows)
    : rows_(rows),
      data_(rows_ * (rows_ + 1) / 2, Scalar{}) {}

template<typename Scalar>
UnboundedLowerTriangularMatrix<Scalar>::UnboundedLowerTriangularMatrix(
    int const rows,
    uninitialized_t)
    : rows_(rows),
      data_(rows_ * (rows_ + 1) / 2) {}

template<typename Scalar>
UnboundedLowerTriangularMatrix<Scalar>::UnboundedLowerTriangularMatrix(
    std::initializer_list<Scalar> data)
    : rows_(static_cast<int>(std::lround((-1 + Sqrt(8 * data.size())) * 0.5))),
      data_(std::move(data)) {
  DCHECK_EQ(data_.size(), rows_ * (rows_ + 1) / 2);
}

template<typename Scalar>
void UnboundedLowerTriangularMatrix<Scalar>::Extend(int const extra_rows) {
  rows_ += extra_rows;
  data_.resize(rows_ * (rows_ + 1) / 2, Scalar{});
}

template<typename Scalar>
void UnboundedLowerTriangularMatrix<Scalar>::Extend(int const extra_rows,
                                                    uninitialized_t) {
  rows_ += extra_rows;
  data_.resize(rows_ * (rows_ + 1) / 2);
}

template<typename Scalar>
void UnboundedLowerTriangularMatrix<Scalar>::Extend(
    std::initializer_list<Scalar> data) {
  std::move(data.begin(), data.end(), std::back_inserter(data_));
  rows_ = static_cast<int>(std::lround((-1 + Sqrt(8 * data_.size())) * 0.5));
  DCHECK_EQ(data_.size(), rows_ * (rows_ + 1) / 2);
}

template<typename Scalar>
void UnboundedLowerTriangularMatrix<Scalar>::EraseToEnd(
    int const begin_row_index) {
  rows_ = begin_row_index;
  data_.erase(data_.begin() + begin_row_index * (begin_row_index + 1) / 2,
              data_.end());
}

template<typename Scalar>
UnboundedUpperTriangularMatrix<Scalar>
UnboundedLowerTriangularMatrix<Scalar>::Transpose() const {
  UnboundedUpperTriangularMatrix<Scalar> u(rows_, uninitialized);
  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j <= i; ++j) {
      u[j][i] = (*this)[i][j];
    }
  }
  return u;
}

template<typename Scalar>
int UnboundedLowerTriangularMatrix<Scalar>::rows() const {
  return rows_;
}

template<typename Scalar>
int UnboundedLowerTriangularMatrix<Scalar>::dimension() const {
  return data_.size();
}

template<typename Scalar>
bool UnboundedLowerTriangularMatrix<Scalar>::operator==(
    UnboundedLowerTriangularMatrix const& right) const {
  return rows_ == right.rows_ && data_ == right.data_;
}

template<typename Scalar>
Scalar* UnboundedLowerTriangularMatrix<Scalar>::operator[](int const index) {
  DCHECK_LT(index, rows_);
  return &data_[index * (index + 1) / 2];
}

template<typename Scalar>
Scalar const* UnboundedLowerTriangularMatrix<Scalar>::operator[](
    int const index) const {
  DCHECK_LT(index, rows_);
  return &data_[index * (index + 1) / 2];
}

template<typename Scalar>
UnboundedUpperTriangularMatrix<Scalar>::UnboundedUpperTriangularMatrix(
    int const columns)
    : columns_(columns),
      data_(columns_ * (columns_ + 1) / 2, Scalar{}) {}

template<typename Scalar>
UnboundedUpperTriangularMatrix<Scalar>::UnboundedUpperTriangularMatrix(
    int const columns,
    uninitialized_t)
    : columns_(columns),
      data_(columns_ * (columns_ + 1) / 2) {}

template<typename Scalar>
UnboundedUpperTriangularMatrix<Scalar>::UnboundedUpperTriangularMatrix(
    std::initializer_list<Scalar> const& data)
    : columns_(
          static_cast<int>(std::lround((-1 + Sqrt(8 * data.size())) * 0.5))),
      data_(Transpose(data,
                      /*current_columns=*/0,
                      /*extra_columns=*/columns_)) {
  DCHECK_EQ(data_.size(), columns_ * (columns_ + 1) / 2);
}

template<typename Scalar>
void UnboundedUpperTriangularMatrix<Scalar>::Extend(int const extra_columns) {
  columns_ += extra_columns;
  data_.resize(columns_ * (columns_ + 1) / 2, Scalar{});
}

template<typename Scalar>
void UnboundedUpperTriangularMatrix<Scalar>::Extend(int const extra_columns,
                                                    uninitialized_t) {
  columns_ += extra_columns;
  data_.resize(columns_ * (columns_ + 1) / 2);
}

template<typename Scalar>
void UnboundedUpperTriangularMatrix<Scalar>::Extend(
    std::initializer_list<Scalar> const& data) {
  int const new_columns = static_cast<int>(
      std::lround((-1 + Sqrt(8 * (data_.size() + data.size()))) * 0.5));
  auto transposed_data = Transpose(data,
                                   /*current_columns=*/columns_,
                                   /*extra_columns=*/new_columns - columns_);
  columns_ = new_columns;
  std::move(transposed_data.begin(),
            transposed_data.end(),
            std::back_inserter(data_));
  DCHECK_EQ(data_.size(), columns_ * (columns_ + 1) / 2);
}

template<typename Scalar>
void UnboundedUpperTriangularMatrix<Scalar>::EraseToEnd(
    int const begin_column_index) {
  columns_ = begin_column_index;
  data_.erase(data_.begin() + begin_column_index * (begin_column_index + 1) / 2,
              data_.end());
}

template<typename Scalar>
UnboundedLowerTriangularMatrix<Scalar>
UnboundedUpperTriangularMatrix<Scalar>::Transpose() const {
  UnboundedLowerTriangularMatrix<Scalar> l(columns_, uninitialized);
  for (int i = 0; i < columns_; ++i) {
    for (int j = i; j < columns_; ++j) {
      l[j][i] = (*this)[i][j];
    }
  }
  return l;
}

template<typename Scalar>
int UnboundedUpperTriangularMatrix<Scalar>::columns() const {
  return columns_;
}

template<typename Scalar>
int UnboundedUpperTriangularMatrix<Scalar>::dimension() const {
  return data_.size();
}

template<typename Scalar>
bool UnboundedUpperTriangularMatrix<Scalar>::operator==(
    UnboundedUpperTriangularMatrix const& right) const {
  return columns_ == right.columns_ && data_ == right.data_;
}

template<typename Scalar>
template<typename Matrix>
Scalar& UnboundedUpperTriangularMatrix<Scalar>::Row<Matrix>::operator[](
    int const column) {
  DCHECK_LT(column, matrix_.columns_);
  return matrix_.data_[column * (column + 1) / 2 + row_];
}

template<typename Scalar>
template<typename Matrix>
Scalar const&
UnboundedUpperTriangularMatrix<Scalar>::Row<Matrix>::operator[](
    int const column) const {
  DCHECK_LT(column, matrix_.columns_);
  return matrix_.data_[column * (column + 1) / 2 + row_];
}

template<typename Scalar>
template<typename Matrix>
UnboundedUpperTriangularMatrix<Scalar>::Row<Matrix>::Row(Matrix& matrix,
                                                         int const row)
    : matrix_(const_cast<std::remove_const_t<Matrix>&>(matrix)),
      row_(row) {}

template<typename Scalar>
auto UnboundedUpperTriangularMatrix<Scalar>::operator[](int const row)
    -> Row<UnboundedUpperTriangularMatrix> {
  return Row<UnboundedUpperTriangularMatrix>{*this, row};
}

template<typename Scalar>
auto UnboundedUpperTriangularMatrix<Scalar>::operator[](int const row) const
    -> Row<UnboundedUpperTriangularMatrix const> {
  return Row<UnboundedUpperTriangularMatrix const>{*this, row};
}

template<typename Scalar>
auto
UnboundedUpperTriangularMatrix<Scalar>::Transpose(
    std::initializer_list<Scalar> const& data,
    int const current_columns,
    int const extra_columns) ->
  std::vector<Scalar, uninitialized_allocator<Scalar>> {
  // |data| is a trapezoidal slice at the end of the matrix.  This is
  // inconvenient to index, so we start by constructing a rectangular array with
  // |extra_columns| columns and |current_columns + extra_columns| rows padded
  // with junk.
  std::vector<Scalar, uninitialized_allocator<Scalar>> padded;
  {
    padded.reserve(2 * data.size());  // An overestimate.
    int row = 0;
    int column = 0;
    for (auto it = data.begin(); it != data.end();) {
      if (row <= current_columns + column) {
        padded.push_back(*it);
        ++it;
      } else {
        padded.emplace_back();
      }
      ++column;
      if (column == extra_columns) {
        column = 0;
        ++row;
      }
    }
  }

  // Scan the padded array by column and append the part above the diagonal to
  // the result.
  std::vector<Scalar, uninitialized_allocator<Scalar>> result;
  result.reserve(data.size());
  int const number_of_rows = current_columns + extra_columns;
  for (int column = 0; column < extra_columns; ++column) {
    for (int row = 0; row < number_of_rows; ++row) {
      if (row <= current_columns + column) {
        result.push_back(padded[row * extra_columns + column]);
      }
    }
  }

  return result;
}

// [Hig02], Algorithm 10.2.
template<typename Scalar>
UnboundedUpperTriangularMatrix<SquareRoot<Scalar>> CholeskyDecomposition(
    UnboundedUpperTriangularMatrix<Scalar> const& A) {
  UnboundedUpperTriangularMatrix<SquareRoot<Scalar>> R(A.columns(),
                                                       uninitialized);
  for (int j = 0; j < A.columns(); ++j) {
    for (int i = 0; i < j; ++i) {
      Scalar Σrₖᵢrₖⱼ{};
      for (int k = 0; k < i; ++k) {
        Σrₖᵢrₖⱼ += R[k][i] * R[k][j];
      }
      R[i][j] = (A[i][j] - Σrₖᵢrₖⱼ) / R[i][i];
    }
    Scalar Σrₖⱼ²{};
    for (int k = 0; k < j; ++k) {
      Σrₖⱼ² += Pow<2>(R[k][j]);
    }
    // This will produce NaNs if the matrix is not positive definite.
    R[j][j] = Sqrt(A[j][j] - Σrₖⱼ²);
  }
  return R;
}

// [KM13], formulæ (10) and (11).
template<typename Scalar>
void ᵗRDRDecomposition(UnboundedUpperTriangularMatrix<Scalar> const& A,
                       UnboundedUpperTriangularMatrix<double>& R,
                       UnboundedVector<Scalar>& D) {
  for (int i = 0; i < A.columns(); ++i) {
    Scalar Σrₖᵢ²dₖ{};
    for (int k = 0; k < i; ++k) {
      Σrₖᵢ²dₖ += Pow<2>(R[k][i]) * D[k];
    }
    D[i] = A[i][i] - Σrₖᵢ²dₖ;
    for (int j = i + 1; j < A.columns(); ++j) {
      Scalar Σrₖᵢrₖⱼdₖ{};
      for (int k = 0; k < i; ++k) {
        Σrₖᵢrₖⱼdₖ += R[k][i] * R[k][j] * D[k];
      }
      R[i][j] = (A[i][j] - Σrₖᵢrₖⱼdₖ) / D[i];
    }
    R[i][i] = 1;
  }
}

// [Hig02], Algorithm 8.1.
template<typename LScalar, typename RScalar>
UnboundedVector<Quotient<RScalar, LScalar>> BackSubstitution(
    UnboundedUpperTriangularMatrix<LScalar> const& U,
    UnboundedVector<RScalar> const& b) {
  UnboundedVector<Quotient<RScalar, LScalar>> x(b.size(), uninitialized);
  int const n = b.size() - 1;
  x[n] = b[n] / U[n][n];
  for (int i = n - 1; i >= 0; --i) {
    auto s = b[i];
    for (int j = i + 1; j <= n; ++j) {
      s -= U[i][j] * x[j];
    }
    x[i] = s / U[i][i];
  }
  return x;
}

// [Hig02] says: "We will not state the analogous algorithm for solving a lower
// triangular system, forward substitution."  So we follow
// https://en.wikipedia.org/wiki/Triangular_matrix#Forward_substitution.
template<typename LScalar, typename RScalar>
UnboundedVector<Quotient<RScalar, LScalar>> ForwardSubstitution(
    UnboundedLowerTriangularMatrix<LScalar> const& L,
    UnboundedVector<RScalar> const& b) {
  UnboundedVector<Quotient<RScalar, LScalar>> x(b.size(), uninitialized);
  x[0] = b[0] / L[0][0];
  for (int i = 1; i < b.size(); ++i) {
    auto s = b[i];
    for (int j = 0; j < i; ++j) {
      s -= L[i][j] * x[j];
    }
    x[i] = s / L[i][i];
  }
  return x;
}

template<typename Scalar>
std::ostream& operator<<(std::ostream& out,
                         UnboundedVector<Scalar> const& vector) {
  std::stringstream s;
  for (int i = 0; i < vector.data_.size(); ++i) {
    s << (i == 0 ? "{" : "") << vector.data_[i]
      << (i == vector.data_.size() - 1 ? "}" : ", ");
  }
  out << s.str();
  return out;
}

template<typename Scalar>
std::ostream& operator<<(std::ostream& out,
                         UnboundedLowerTriangularMatrix<Scalar> const& matrix) {
  out << "rows: " << matrix.rows() << "\n";
  for (int i = 0; i < matrix.rows(); ++i) {
    out << "{";
    for (int j = 0; j <= i; ++j) {
      out << matrix[i][j];
      if (j < i) {
        out << ", ";
      }
    }
    out << "}\n";
  }
  return out;
}

template<typename Scalar>
std::ostream& operator<<(std::ostream& out,
                         UnboundedUpperTriangularMatrix<Scalar> const& matrix) {
  out << "columns: " << matrix.columns_ << "\n";
  for (int i = 0; i < matrix.columns(); ++i) {
    out << "{";
    for (int j = i; j < matrix.columns(); ++j) {
      if (j > i) {
        out << ", ";
      }
      out << matrix[i][j];
    }
    out << "}\n";
  }
  return out;
}

}  // namespace internal_unbounded_arrays
}  // namespace numerics
}  // namespace principia
