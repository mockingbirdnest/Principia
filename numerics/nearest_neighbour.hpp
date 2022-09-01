#pragma once

#include <vector>

namespace principia {
namespace numerics {
namespace internal_nearest_neighbour {

template<typename Value_>
class PrincipalComponentPartitioningTree {
 public:
  using Value = Value_;

  PrincipalComponentPartitioningTree(std::vector<Value> const& values);

  void Add(Value const& value);

  Value const& FindNearest(Value const& value) const;

 private:
  Value centroid_;
};

}  // namespace internal_nearest_neighbour

using internal_nearest_neighbour::PrincipalComponentPartitioningTree;

}  // namespace numerics
}  // namespace principia

#include "numerics/nearest_neighbour_body.hpp"