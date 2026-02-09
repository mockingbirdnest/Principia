#pragma once

#include <array>
#include <utility>
#include <vector>

#include "base/algebra.hpp"
#include "base/tags.hpp"
#include "numerics/matrix_views.hpp"
#include "numerics/transposed_view.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace numerics {
namespace _fixed_arrays {
namespace internal {

using namespace principia::base::_algebra;
using namespace principia::base::_tags;
using namespace principia::numerics::_matrix_views;
using namespace principia::numerics::_transposed_view;
using namespace principia::quantities::_si;

//TODO(phl)Conversions across use_heap?
template<typename Scalar_, std::int64_t rows_, std::int64_t columns_,
         bool use_heap = false>
class FixedMatrix;
template<typename Scalar_, std::int64_t rows_, bool use_heap = false>
class FixedUpperTriangularMatrix;

template<typename Scalar_, std::int64_t size_, bool use_heap = false>
class FixedVector final {
 public:
  using Scalar = Scalar_;
  static constexpr std::int64_t dimension = size_;

  constexpr FixedVector();
  explicit FixedVector(uninitialized_t);

  constexpr FixedVector(
      std::array<Scalar, size_> const& data);  // NOLINT(runtime/explicit)
  constexpr FixedVector(
      std::array<Scalar, size_>&& data);  // NOLINT(runtime/explicit)

  template<bool uh>
  constexpr FixedVector(FixedVector<Scalar, size_, uh> const& other);
  template<bool uh>
  constexpr FixedVector(FixedVector<Scalar, size_, uh>&& other);

  // Constructs a fixed vector by copying data from the view.  Note that the
  // result is fixed even if the matrix being viewed is an UnboundedMatrix.
  template<typename T>
    requires std::same_as<typename T::Scalar, Scalar_>
  constexpr explicit FixedVector(ColumnView<T> const& view);

  // Convertible to an array.
  explicit constexpr operator std::array<Scalar, size_>() const;

  constexpr Scalar& operator[](std::int64_t index);
  constexpr Scalar const& operator[](std::int64_t index) const;

  template<bool uh>
  constexpr FixedVector& operator=(FixedVector<Scalar, size_, uh> const& other);
  template<bool uh>
  constexpr FixedVector& operator=(FixedVector<Scalar, size_, uh>&& other);

  constexpr FixedVector& operator=(Scalar const (&right)[size_]);

  constexpr FixedVector& operator+=(FixedVector const& right);
  constexpr FixedVector& operator-=(FixedVector const& right);
  constexpr FixedVector& operator*=(double right);
  constexpr FixedVector& operator/=(double right);

  Scalar Norm() const;
  Square<Scalar> Norm²() const;

  template<bool uh = use_heap>
  FixedVector<double, size_, uh> Normalize() const;

  static constexpr std::int64_t size() { return size_; }

  typename std::array<Scalar, size_>::const_iterator begin() const;
  typename std::array<Scalar, size_>::const_iterator end() const;

  template<typename H>
  friend H AbslHashValue(H h, FixedVector const& vector) {
    for (std::int64_t index = 0; index < size_; ++index) {
      h = H::combine(std::move(h), vector.data_[index] / si::Unit<Scalar>);
    }
    return h;
  }

 private:
  using Data = std::conditional_t<use_heap,
                                  std::unique_ptr<std::array<Scalar, size_>>,
                                  std::array<Scalar, size_>>;

  std::array<Scalar, size_> const& data() const;
  std::array<Scalar, size_>& data();

  Data data_;

  template<typename S, std::int64_t s, bool luh, bool ruh>
  friend bool operator==(FixedVector<S, s, luh> const& left,
                         FixedVector<S, s, ruh> const& right);
  template<typename L, typename R, std::int64_t s, bool ruh>
  friend constexpr Product<L, R> operator*(
      L* const left,
      FixedVector<R, s, ruh> const& right);
  template<typename L, typename R, std::int64_t s, bool luh, bool ruh>
  friend constexpr Product<L, R> operator*(
      TransposedView<FixedVector<L, s, luh>> const& left,
      FixedVector<R, s, ruh> const& right);
  template<typename L, typename R, std::int64_t r, std::int64_t c,
           bool uh, bool luh, bool ruh>
  friend constexpr FixedVector<Product<L, R>, r, uh> operator*(
      FixedMatrix<L, r, c, luh> const& left,
      FixedVector<R, c, ruh> const& right);

  template<typename S, std::int64_t s, bool uh>
  friend class FixedVector;
};

template<typename Scalar_, std::int64_t rows_, std::int64_t columns_,
         bool use_heap>
class FixedMatrix final {
  static constexpr std::int64_t size_ = rows_ * columns_;

 public:
  using Scalar = Scalar_;
  static constexpr std::int64_t rows() { return rows_; }
  static constexpr std::int64_t columns() { return columns_; }

  constexpr FixedMatrix();
  explicit FixedMatrix(uninitialized_t);

  // The `data` must be in row-major format.
  constexpr FixedMatrix(
      std::array<Scalar, size_> const& data);  // NOLINT(runtime/explicit)
  constexpr FixedMatrix(
      std::array<Scalar, size_>&& data);  // NOLINT(runtime/explicit)

  template<bool uh>
  constexpr FixedMatrix(FixedMatrix<Scalar, rows_, columns_, uh> const& other);
  template<bool uh>
  constexpr FixedMatrix(FixedMatrix<Scalar, rows_, columns_, uh>&& other);

  template<bool uh>
  constexpr explicit FixedMatrix(
      TransposedView<FixedMatrix<Scalar, columns_, rows_, uh>> const& view);

  // For  0 < i < rows and 0 < j < columns, the entry a_ij is accessed as
  // `a(i, j)`.  if i and j do not satisfy these conditions, the expression
  // `a(i, j)` implies undefined behaviour.
  constexpr Scalar& operator()(std::int64_t row, std::int64_t column);
  constexpr Scalar const& operator()(std::int64_t row,
                                     std::int64_t column) const;

  template<bool uh>
  constexpr FixedMatrix& operator=(
      FixedMatrix<Scalar, rows_, columns_, uh> const& other);
  template<bool uh>
  constexpr FixedMatrix& operator=(
      FixedMatrix<Scalar, rows_, columns_, uh>&& other);

  constexpr FixedMatrix& operator=(Scalar const (&right)[size_]);

  constexpr FixedMatrix& operator+=(FixedMatrix const& right);
  constexpr FixedMatrix& operator-=(FixedMatrix const& right);
  constexpr FixedMatrix& operator*=(double right);
  constexpr FixedMatrix& operator/=(double right);

  constexpr FixedMatrix& operator*=(
      FixedMatrix<double, rows_, columns_, use_heap> const& right)
    requires(rows_ == columns_);

  template<std::int64_t r>
  Scalar const* row() const;

  Scalar FrobeniusNorm() const;

  // Applies the matrix as a bilinear form.  Present for compatibility with
  // `SymmetricBilinearForm`.  Prefer to use `TransposedView` and `operator*`.
  template<typename LScalar, typename RScalar, bool luh, bool ruh>
  Product<Scalar, Product<LScalar, RScalar>>
      operator()(FixedVector<LScalar, columns_, luh> const& left,
                 FixedVector<RScalar, rows_, ruh> const& right) const;

  static FixedMatrix Identity()
    requires(std::is_arithmetic_v<Scalar_> && rows_ == columns_);

 private:
  using Data = std::conditional_t<use_heap,
                                  std::unique_ptr<std::array<Scalar, size_>>,
                                  std::array<Scalar, size_>>;

  std::array<Scalar, size_> const& data() const;
  std::array<Scalar, size_>& data();

  Data data_;

  template<typename S, std::int64_t r, std::int64_t c, bool luh, bool ruh>
  friend bool operator==(FixedMatrix<S, r, c, luh> const& left,
                         FixedMatrix<S, r, c, ruh> const& right);
  template<typename L, typename R, std::int64_t r, std::int64_t c,
           bool uh, bool luh, bool ruh>
  friend constexpr FixedVector<Product<L, R>, r, uh> operator*(
      FixedMatrix<L, r, c, luh> const& left,
      FixedVector<R, c, ruh> const& right);

  template<typename S, std::int64_t r, std::int64_t c, bool uh>
  friend class FixedMatrix;
};

template<typename Scalar_, std::int64_t rows_, bool use_heap = false>
class FixedStrictlyLowerTriangularMatrix final {
  static constexpr std::int64_t size_ = rows_ * (rows_ - 1) / 2;

 public:
  using Scalar = Scalar_;
  static constexpr std::int64_t rows() { return rows_; }
  static constexpr std::int64_t columns() { return rows_; }

  constexpr FixedStrictlyLowerTriangularMatrix();
  explicit FixedStrictlyLowerTriangularMatrix(uninitialized_t);

  // The `data` must be in row-major format.
  constexpr FixedStrictlyLowerTriangularMatrix(
      std::array<Scalar, size_> const& data);

  template<bool uh>
  constexpr FixedStrictlyLowerTriangularMatrix(
      FixedStrictlyLowerTriangularMatrix<Scalar, rows_, uh> const& other);
  template<bool uh>
  constexpr FixedStrictlyLowerTriangularMatrix(
      FixedStrictlyLowerTriangularMatrix<Scalar, rows_, uh>&& other);

  template<bool uh = use_heap>
  explicit constexpr operator FixedMatrix<Scalar, rows_, rows_, uh>() const;

  // For  0 ≤ j < i < rows, the entry a_ij is accessed as `a(i, j)`.
  // if i and j do not satisfy these conditions, the expression `a(i, j)`
  // implies undefined behaviour.
  constexpr Scalar& operator()(std::int64_t row, std::int64_t column);
  constexpr Scalar const& operator()(std::int64_t row,
                                     std::int64_t column) const;

  template<bool uh>
  constexpr FixedStrictlyLowerTriangularMatrix& operator=(
      FixedStrictlyLowerTriangularMatrix<Scalar, rows_, uh> const& other);
  template<bool uh>
  constexpr FixedStrictlyLowerTriangularMatrix& operator=(
      FixedStrictlyLowerTriangularMatrix<Scalar, rows_, uh>&& other);

  constexpr FixedStrictlyLowerTriangularMatrix& operator=(
      Scalar const (&right)[size_]);

  template<std::int64_t r>
  Scalar const* row() const;

 private:
  using Data = std::conditional_t<use_heap,
                                  std::unique_ptr<std::array<Scalar, size_>>,
                                  std::array<Scalar, size_>>;

  std::array<Scalar, size_> const& data() const;
  std::array<Scalar, size_>& data();

  Data data_;

  template<typename S, std::int64_t r,  bool luh, bool ruh>
  friend bool operator==(
      FixedStrictlyLowerTriangularMatrix<S, r, luh> const& left,
      FixedStrictlyLowerTriangularMatrix<S, r, ruh> const& right);

  template<typename S, std::int64_t r, bool uh>
  friend class FixedStrictlyLowerTriangularMatrix;
};

template<typename Scalar_, std::int64_t rows_, bool use_heap = false>
class FixedLowerTriangularMatrix final {
  static constexpr std::int64_t size_ = rows_ * (rows_ + 1) / 2;

 public:
  using Scalar = Scalar_;
  static constexpr std::int64_t rows() { return rows_; }
  static constexpr std::int64_t columns() { return rows_; }

  constexpr FixedLowerTriangularMatrix();
  explicit FixedLowerTriangularMatrix(uninitialized_t);

  // The `data` must be in row-major format.
  constexpr FixedLowerTriangularMatrix(
      std::array<Scalar, size_> const& data);

  template<bool uh>
  constexpr FixedLowerTriangularMatrix(
      FixedLowerTriangularMatrix<Scalar, rows_, uh> const& other);
  template<bool uh>
  constexpr FixedLowerTriangularMatrix(
      FixedLowerTriangularMatrix<Scalar, rows_, uh>&& other);

  template<bool uh>
  explicit FixedLowerTriangularMatrix(
      TransposedView<FixedUpperTriangularMatrix<Scalar, rows_, uh>> const&
          view);

  template<bool uh = use_heap>
  explicit constexpr operator FixedMatrix<Scalar, rows_, rows_, uh>() const;

  // For  0 ≤ j ≤ i < rows, the entry a_ij is accessed as `a(i, j)`.
  // if i and j do not satisfy these conditions, the expression `a(i, j)`
  // implies undefined behaviour.
  constexpr Scalar& operator()(std::int64_t row, std::int64_t column);
  constexpr Scalar const& operator()(std::int64_t row,
                                     std::int64_t column) const;

  template<bool uh>
  constexpr FixedLowerTriangularMatrix& operator=(
      FixedLowerTriangularMatrix<Scalar, rows_, uh> const& other);
  template<bool uh>
  constexpr FixedLowerTriangularMatrix& operator=(
      FixedLowerTriangularMatrix<Scalar, rows_, uh>&& other);

  constexpr FixedLowerTriangularMatrix& operator=(
      Scalar const (&right)[size_]);

 private:
  using Data = std::conditional_t<use_heap,
                                  std::unique_ptr<std::array<Scalar, size_>>,
                                  std::array<Scalar, size_>>;

  std::array<Scalar, size_> const& data() const;
  std::array<Scalar, size_>& data();

  Data data_;

  template<typename S, std::int64_t r,  bool luh, bool ruh>
  friend bool operator==(FixedLowerTriangularMatrix<S, r, luh> const& left,
                         FixedLowerTriangularMatrix<S, r, ruh> const& right);

  template<typename S, std::int64_t r, bool uh>
  friend class FixedLowerTriangularMatrix;
};

template<typename Scalar_, std::int64_t columns_, bool use_heap = false>
class FixedStrictlyUpperTriangularMatrix final {
  static constexpr std::int64_t size_ = columns_ * (columns_ - 1) / 2;

 public:
  using Scalar = Scalar_;
  static constexpr std::int64_t rows() { return columns_; }
  static constexpr std::int64_t columns() { return columns_; }

  constexpr FixedStrictlyUpperTriangularMatrix();
  explicit FixedStrictlyUpperTriangularMatrix(uninitialized_t);

  // The `data` must be in row-major format.
  constexpr FixedStrictlyUpperTriangularMatrix(
      std::array<Scalar, size_> const& data);

  template<bool uh>
  constexpr FixedStrictlyUpperTriangularMatrix(
      FixedStrictlyUpperTriangularMatrix<Scalar, columns_, uh> const& other);
  template<bool uh>
  constexpr FixedStrictlyUpperTriangularMatrix(
      FixedStrictlyUpperTriangularMatrix<Scalar, columns_, uh>&& other);

  template<bool uh>
  explicit FixedStrictlyUpperTriangularMatrix(
      TransposedView<
          FixedStrictlyLowerTriangularMatrix<Scalar, columns_, uh>> const&
          view);

  template<bool uh = use_heap>
  explicit constexpr
  operator FixedMatrix<Scalar, columns_, columns_, uh>() const;

  // For  0 ≤ i < j < columns, the entry a_ij is accessed as `a(i, j)`.
  // if i and j do not satisfy these conditions, the expression `a(i, j)`
  // implies undefined behaviour.
  constexpr Scalar& operator()(std::int64_t row, std::int64_t column);
  constexpr Scalar const& operator()(std::int64_t row,
                                     std::int64_t column) const;

  template<bool uh>
  constexpr FixedStrictlyUpperTriangularMatrix& operator=(
      FixedStrictlyUpperTriangularMatrix<Scalar, columns_, uh> const& other);
  template<bool uh>
  constexpr FixedStrictlyUpperTriangularMatrix& operator=(
      FixedStrictlyUpperTriangularMatrix<Scalar, columns_, uh>&& other);

  constexpr FixedStrictlyUpperTriangularMatrix& operator=(
      Scalar const (&right)[size_]);

 private:
  using Data = std::conditional_t<use_heap,
                                  std::unique_ptr<std::array<Scalar, size_>>,
                                  std::array<Scalar, size_>>;

  // For ease of writing matrices in tests, the input data is received in row-
  // major format.  This transposes a trapezoidal slice to make it column-major.
  static std::array<Scalar, size_> Transpose(
      std::array<Scalar, size_> const& data);

  std::array<Scalar, size_> const& data() const;
  std::array<Scalar, size_>& data();

  Data data_;

  template<typename S, std::int64_t c,  bool luh, bool ruh>
  friend bool operator==(
      FixedStrictlyUpperTriangularMatrix<S, c, luh> const& left,
      FixedStrictlyUpperTriangularMatrix<S, c, ruh> const& right);

  template<typename S, std::int64_t c, bool uh>
  friend class FixedStrictlyUpperTriangularMatrix;
};

template<typename Scalar_, std::int64_t columns_, bool use_heap>
class FixedUpperTriangularMatrix final {
  static constexpr std::int64_t size_ = columns_ * (columns_ + 1) / 2;

 public:
  using Scalar = Scalar_;
  static constexpr std::int64_t rows() { return columns_; }
  static constexpr std::int64_t columns() { return columns_; }

  constexpr FixedUpperTriangularMatrix();
  explicit FixedUpperTriangularMatrix(uninitialized_t);

  // The `data` must be in row-major format.
  constexpr FixedUpperTriangularMatrix(
      std::array<Scalar, size_> const& data);

  template<bool uh>
  constexpr FixedUpperTriangularMatrix(
      FixedUpperTriangularMatrix<Scalar, columns_, uh> const& other);
  template<bool uh>
  constexpr FixedUpperTriangularMatrix(
      FixedUpperTriangularMatrix<Scalar, columns_, uh>&& other);

  template<bool uh>
  explicit FixedUpperTriangularMatrix(
      TransposedView<
          FixedLowerTriangularMatrix<Scalar, columns_, uh>> const& view);

  template<bool uh = use_heap>
  explicit constexpr
  operator FixedMatrix<Scalar, columns_, columns_, uh>() const;

  // For  0 ≤ i ≤ j < columns, the entry a_ij is accessed as `a(i, j)`.
  // if i and j do not satisfy these conditions, the expression `a(i, j)`
  // implies undefined behaviour.
  constexpr Scalar& operator()(std::int64_t row, std::int64_t column);
  constexpr Scalar const& operator()(std::int64_t row,
                                     std::int64_t column) const;

  template<bool uh>
  constexpr FixedUpperTriangularMatrix& operator=(
      FixedUpperTriangularMatrix<Scalar, columns_, uh> const& other);
  template<bool uh>
  constexpr FixedUpperTriangularMatrix& operator=(
      FixedUpperTriangularMatrix<Scalar, columns_, uh>&& other);

  constexpr FixedUpperTriangularMatrix& operator=(
      Scalar const (&right)[size_]);

 private:
  using Data = std::conditional_t<use_heap,
                                  std::unique_ptr<std::array<Scalar, size_>>,
                                  std::array<Scalar, size_>>;

  // For ease of writing matrices in tests, the input data is received in row-
  // major format.  This transposes a trapezoidal slice to make it column-major.
  static std::array<Scalar, size_> Transpose(
      std::array<Scalar, size_> const& data);

  std::array<Scalar, size_> const& data() const;
  std::array<Scalar, size_>& data();

  Data data_;

  template<typename S, std::int64_t c,  bool luh, bool ruh>
  friend bool operator==(FixedUpperTriangularMatrix<S, c, luh> const& left,
                         FixedUpperTriangularMatrix<S, c, ruh> const& right);

  template<typename S, std::int64_t c, bool uh>
  friend class FixedUpperTriangularMatrix;
};

// Prefer using the operator* that takes a TransposedView.
template<typename LScalar, typename RScalar, std::int64_t size,
         bool luh, bool ruh>
constexpr Product<LScalar, RScalar> InnerProduct(
    FixedVector<LScalar, size, luh> const& left,
    FixedVector<RScalar, size, ruh> const& right);

template<typename Scalar, std::int64_t size, bool uh, bool vuh>
constexpr FixedVector<double, size, uh> Normalize(
    FixedVector<Scalar, size, vuh> const& vector);

template<typename LScalar, typename RScalar, std::int64_t size,
         bool uh, bool luh, bool ruh>
constexpr FixedMatrix<Product<LScalar, RScalar>, size, size, uh>
SymmetricProduct(FixedVector<LScalar, size, luh> const& left,
                 FixedVector<RScalar, size, ruh> const& right);

template<typename Scalar, std::int64_t size, bool uh, bool vuh>
constexpr FixedMatrix<Square<Scalar>, size, size, uh> SymmetricSquare(
    FixedVector<Scalar, size, vuh> const& vector);

// Equality.

template<typename Scalar, std::int64_t size,
         bool luh, bool ruh>
bool operator==(FixedVector<Scalar, size, luh> const& left,
                FixedVector<Scalar, size, ruh> const& right);

template<typename Scalar, std::int64_t rows, std::int64_t columns,
         bool luh, bool ruh>
bool operator==(FixedMatrix<Scalar, rows, columns, luh> const& left,
                FixedMatrix<Scalar, rows, columns, ruh> const& right);

template<typename Scalar, std::int64_t rows,
         bool luh, bool ruh>
bool operator==(
    FixedStrictlyLowerTriangularMatrix<Scalar, rows, luh> const& left,
    FixedStrictlyLowerTriangularMatrix<Scalar, rows, ruh> const& right);

template<typename Scalar, std::int64_t rows,
         bool luh, bool ruh>
bool operator==(FixedLowerTriangularMatrix<Scalar, rows, luh> const& left,
                FixedLowerTriangularMatrix<Scalar, rows, ruh> const& right);

template<typename Scalar, std::int64_t columns,
         bool luh, bool ruh>
bool operator==(
    FixedStrictlyUpperTriangularMatrix<Scalar, columns, luh> const& left,
    FixedStrictlyUpperTriangularMatrix<Scalar, columns, ruh> const& right);

template<typename Scalar, std::int64_t columns,
         bool luh, bool ruh>
bool operator==(FixedUpperTriangularMatrix<Scalar, columns, luh> const& left,
                FixedUpperTriangularMatrix<Scalar, columns, ruh> const& right);

template<typename Scalar, std::int64_t size,
         bool luh, bool ruh>
bool operator!=(FixedVector<Scalar, size, luh> const& left,
                FixedVector<Scalar, size, ruh> const& right);

template<typename Scalar, std::int64_t rows, std::int64_t columns,
         bool luh, bool ruh>
bool operator!=(FixedMatrix<Scalar, rows, columns, luh> const& left,
                FixedMatrix<Scalar, rows, columns, ruh> const& right);

template<typename Scalar, std::int64_t rows,
         bool luh, bool ruh>
bool operator!=(
    FixedStrictlyLowerTriangularMatrix<Scalar, rows, luh> const& left,
    FixedStrictlyLowerTriangularMatrix<Scalar, rows, ruh> const& right);

template<typename Scalar, std::int64_t rows,
         bool luh, bool ruh>
bool operator!=(FixedLowerTriangularMatrix<Scalar, rows, luh> const& left,
                FixedLowerTriangularMatrix<Scalar, rows, ruh> const& right);

template<typename Scalar, std::int64_t columns,
         bool luh, bool ruh>
bool operator!=(
    FixedStrictlyUpperTriangularMatrix<Scalar, columns, luh> const& left,
    FixedStrictlyUpperTriangularMatrix<Scalar, columns, ruh> const& right);

template<typename Scalar, std::int64_t columns,
         bool luh, bool ruh>
bool operator!=(FixedUpperTriangularMatrix<Scalar, columns, luh> const& left,
                FixedUpperTriangularMatrix<Scalar, columns, ruh> const& right);

// Additive groups.

template<typename Scalar, std::int64_t size,
         bool uh = false, bool ruh>
constexpr FixedVector<Scalar, size, uh> operator+(
    FixedVector<Scalar, size, ruh> const& right);

template<typename Scalar, std::int64_t rows, std::int64_t columns,
         bool uh = false, bool ruh>
constexpr FixedMatrix<Scalar, rows, columns, uh> operator+(
    FixedMatrix<Scalar, rows, columns, ruh> const& right);

template<typename Scalar, std::int64_t size,
         bool uh = false, bool ruh>
constexpr FixedVector<Scalar, size, uh> operator-(
    FixedVector<Scalar, size, ruh> const& right);

template<typename Scalar, std::int64_t rows, std::int64_t columns,
         bool uh = false, bool ruh>
constexpr FixedMatrix<Scalar, rows, columns, uh> operator-(
    FixedMatrix<Scalar, rows, columns, ruh> const& right);

template<typename LScalar, typename RScalar, std::int64_t size,
         bool uh = false, bool luh, bool ruh>
constexpr FixedVector<Sum<LScalar, RScalar>, size, uh> operator+(
    FixedVector<LScalar, size, luh> const& left,
    FixedVector<RScalar, size, ruh> const& right);

template<typename LScalar, typename RScalar,
         std::int64_t rows, std::int64_t columns,
         bool uh = false, bool luh, bool ruh>
constexpr FixedMatrix<Sum<LScalar, RScalar>, rows, columns, uh> operator+(
    FixedMatrix<LScalar, rows, columns, luh> const& left,
    FixedMatrix<RScalar, rows, columns, ruh> const& right);

template<typename LScalar, typename RScalar, std::int64_t size,
         bool uh = false, bool luh, bool ruh>
constexpr FixedVector<Difference<LScalar, RScalar>, size, uh> operator-(
    FixedVector<LScalar, size, luh> const& left,
    FixedVector<RScalar, size, ruh> const& right);

template<typename LScalar, typename RScalar,
         std::int64_t rows, std::int64_t columns,
         bool uh = false, bool luh, bool ruh>
constexpr FixedMatrix<Difference<LScalar, RScalar>, rows, columns, uh>
operator-(FixedMatrix<LScalar, rows, columns, luh> const& left,
          FixedMatrix<RScalar, rows, columns, ruh> const& right);

// Vector spaces.

template<typename LScalar, typename RScalar, std::int64_t size,
         bool uh = false, bool ruh>
constexpr FixedVector<Product<LScalar, RScalar>, size, uh> operator*(
    LScalar const& left,
    FixedVector<RScalar, size, ruh> const& right);

template<typename LScalar, typename RScalar, std::int64_t size,
         bool uh = false, bool luh>
constexpr FixedVector<Product<LScalar, RScalar>, size, uh> operator*(
    FixedVector<LScalar, size, luh> const& left,
    RScalar const& right);

template<typename LScalar, typename RScalar,
         std::int64_t rows, std::int64_t columns,
         bool uh = false, bool ruh>
constexpr FixedMatrix<Product<LScalar, RScalar>, rows, columns, uh>
operator*(LScalar const& left,
          FixedMatrix<RScalar, rows, columns, ruh> const& right);

template<typename LScalar, typename RScalar,
         std::int64_t rows, std::int64_t columns,
         bool uh = false, bool luh>
constexpr FixedMatrix<Product<LScalar, RScalar>, rows, columns, uh>
operator*(FixedMatrix<LScalar, rows, columns, luh> const& left,
          RScalar const& right);

template<typename LScalar, typename RScalar, std::int64_t size,
         bool uh = false, bool luh>
constexpr FixedVector<Quotient<LScalar, RScalar>, size, uh> operator/(
    FixedVector<LScalar, size, luh> const& left,
    RScalar const& right);

template<typename LScalar, typename RScalar,
         std::int64_t rows, std::int64_t columns,
         bool uh = false, bool luh>
constexpr FixedMatrix<Quotient<LScalar, RScalar>, rows, columns, uh>
operator/(FixedMatrix<LScalar, rows, columns, luh> const& left,
          RScalar const& right);

// Hilbert space and algebra.

// TODO(phl): We should have a RowView.
template<typename LScalar, typename RScalar, std::int64_t size, bool ruh>
constexpr Product<LScalar, RScalar> operator*(
    LScalar* const left,
    FixedVector<RScalar, size, ruh> const& right);

template<typename LScalar, typename RScalar, std::int64_t size,
         bool luh, bool ruh>
constexpr Product<LScalar, RScalar> operator*(
    TransposedView<FixedVector<LScalar, size, luh>> const& left,
    FixedVector<RScalar, size, ruh> const& right);

template<typename LScalar, typename RScalar,
         std::int64_t lsize, std::int64_t rsize,
         bool uh = false, bool luh, bool ruh>
constexpr FixedMatrix<Product<LScalar, RScalar>, lsize, rsize, uh> operator*(
    FixedVector<LScalar, lsize, luh> const& left,
    TransposedView<FixedVector<RScalar, rsize, ruh>> const& right);

template<typename LScalar, typename RScalar,
         std::int64_t rows, std::int64_t dimension, std::int64_t columns,
         bool uh = false, bool luh, bool ruh>
constexpr FixedMatrix<Product<LScalar, RScalar>, rows, columns, uh>
operator*(FixedMatrix<LScalar, rows, dimension, luh> const& left,
          FixedMatrix<RScalar, dimension, columns, ruh> const& right);

template<typename LScalar, typename RScalar,
         std::int64_t rows, std::int64_t columns,
         bool uh = false, bool luh, bool ruh>
constexpr FixedVector<Product<LScalar, RScalar>, rows, uh> operator*(
    FixedMatrix<LScalar, rows, columns, luh> const& left,
    FixedVector<RScalar, columns, ruh> const& right);

// Use this operator to multiply a row vector with a matrix.  We don't have an
// operator returning a TransposedView as that would cause dangling references.
template<typename LScalar, typename RScalar,
         std::int64_t rows, std::int64_t columns,
         bool uh = false, bool luh, bool ruh>
constexpr FixedVector<Product<LScalar, RScalar>, columns, uh> operator*(
    TransposedView<FixedMatrix<LScalar, rows, columns, luh>> const& left,
    FixedVector<RScalar, rows, ruh> const& right);

// Ouput.

template<typename Scalar, std::int64_t size, bool uh>
std::ostream& operator<<(std::ostream& out,
                         FixedVector<Scalar, size, uh> const& vector);

template<typename Scalar, std::int64_t rows, std::int64_t columns, bool uh>
std::ostream& operator<<(std::ostream& out,
                         FixedMatrix<Scalar, rows, columns, uh> const& matrix);

template<typename Scalar, std::int64_t rows, bool uh>
std::ostream& operator<<(
    std::ostream& out,
    FixedStrictlyLowerTriangularMatrix<Scalar, rows, uh> const& matrix);

template<typename Scalar, std::int64_t rows, bool uh>
std::ostream& operator<<(
    std::ostream& out,
    FixedLowerTriangularMatrix<Scalar, rows, uh> const& matrix);

template<typename Scalar, std::int64_t columns, bool uh>
std::ostream& operator<<(
    std::ostream& out,
    FixedUpperTriangularMatrix<Scalar, columns, uh> const& matrix);

}  // namespace internal

using internal::FixedLowerTriangularMatrix;
using internal::FixedMatrix;
using internal::FixedStrictlyLowerTriangularMatrix;
using internal::FixedStrictlyUpperTriangularMatrix;
using internal::FixedUpperTriangularMatrix;
using internal::FixedVector;

}  // namespace _fixed_arrays
}  // namespace numerics
}  // namespace principia

#include "numerics/fixed_arrays_body.hpp"
