
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
  ::new ((void*)p) U(std::forward<Args>(args)...);
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

}  // namespace internal_unbounded_arrays
}  // namespace numerics
}  // namespace principia
