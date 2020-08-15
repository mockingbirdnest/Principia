
#pragma once

#include <vector>

#include "base/tags.hpp"

namespace principia {
namespace numerics {
namespace internal_unbounded_arrays {

using base::uninitialized_t;

template<typename Scalar>
class UnboundedLowerTriangularMatrix final {
 public:
  UnboundedLowerTriangularMatrix(int rows);
  UnboundedLowerTriangularMatrix(int rows, uninitialized_t);

  void Extend(int rows);
  void Extend(int rows, uninitialized_t);

  int rows() const;
  int dimension() const;

  bool operator==(UnboundedLowerTriangularMatrix const& right) const;

  // For  0 < j <= i < rows, the entry a_ij is accessed as |a[i][j]|.
  // if i and j do not satisfy these conditions, the expression |a[i][j]| is
  // erroneous.
  Scalar* operator[](int index);
  Scalar const* operator[](int index) const;

 private:
  int rows_;
  std::vector<Scalar> data_;
};

}  // namespace internal_unbounded_arrays

using internal_unbounded_arrays::UnboundedLowerTriangularMatrix;

}  // namespace numerics
}  // namespace principia

#include "numerics/unbounded_arrays_body.hpp"
