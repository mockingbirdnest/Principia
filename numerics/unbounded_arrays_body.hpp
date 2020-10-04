
#pragma once

#include "numerics/unbounded_arrays.hpp"

#include <cmath>

#include "base/macros.hpp"
#include "quantities/elementary_functions.hpp"
#include "unbounded_arrays.hpp"

namespace principia {
namespace numerics {
namespace internal_unbounded_arrays {

using base::uninitialized;
using quantities::Sqrt;

template<class T>
template<class U, class... Args>
inline void uninitialized_allocator<T>::construct(U* const p, Args&&... args) {
#if PRINCIPIA_COMPILER_CLANG
  ::new ((void*)p) U(std::forward<Args>(args)...);  // NOLINT
#endif
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
template<typename Reference>
Scalar& UnboundedUpperTriangularMatrix<Scalar>::Row<Reference>::operator[](
    int const column) {
  DCHECK_LT(column, reference_.columns_);
  return reference_.data_[column * (column + 1) / 2 + row_];
}

template<typename Scalar>
template<typename Reference>
Scalar const&
UnboundedUpperTriangularMatrix<Scalar>::Row<Reference>::operator[](
    int const column) const {
  DCHECK_LT(column, reference_.columns_);
  return reference_.data_[column * (column + 1) / 2 + row_];
}

template<typename Scalar>
template<typename Reference>
UnboundedUpperTriangularMatrix<Scalar>::Row<Reference>::Row(Reference reference,
                                                            int const row)
    : reference_(reference),
      row_(row) {}

template<typename Scalar>
auto UnboundedUpperTriangularMatrix<Scalar>::operator[](int const row)
    -> Row<UnboundedUpperTriangularMatrix&> {
  return Row<UnboundedUpperTriangularMatrix&>{*this, row};
}

template<typename Scalar>
auto UnboundedUpperTriangularMatrix<Scalar>::operator[](int const row) const
    -> Row<UnboundedUpperTriangularMatrix const&> {
  return Row<UnboundedUpperTriangularMatrix const&>{*this, row};
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
  std::stringstream s;
  // TODO(phl): Triangular printout.
  s << "rows: " << matrix.rows_ << " ";
  for (int i = 0; i < matrix.data_.size(); ++i) {
    s << (i == 0 ? "{" : "") << matrix.data_[i]
      << (i == matrix.data_.size() - 1 ? "}" : ", ");
  }
  out << s.str();
  return out;
}

template<typename Scalar>
std::ostream& operator<<(std::ostream& out,
                         UnboundedUpperTriangularMatrix<Scalar> const& matrix) {
  std::stringstream s;
  // TODO(phl): Triangular printout.
  s << "columns: " << matrix.columns_ << " ";
  for (int i = 0; i < matrix.data_.size(); ++i) {
    s << (i == 0 ? "{" : "") << matrix.data_[i]
      << (i == matrix.data_.size() - 1 ? "}" : ", ");
  }
  out << s.str();
  return out;
}

}  // namespace internal_unbounded_arrays
}  // namespace numerics
}  // namespace principia
