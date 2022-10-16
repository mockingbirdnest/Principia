#pragma once

#include "numerics/global_optimization.hpp"

namespace principia {
namespace numerics {
namespace internal_global_optimization {

template<typename Scalar, typename Frame>
MultiLevelSingleLinkage<Scalar, Frame>::MultiLevelSingleLinkage(
    Box const& box,
    Field<Scalar, Frame> const& f,
    Field<Gradient<Scalar, Frame>, Frame> const& grad_f,
    int64_t const values_per_round) {}

}  // namespace internal_global_optimization
}  // namespace numerics
}  // namespace principia
