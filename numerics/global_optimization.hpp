#pragma once

#include <functional>

#include "geometry/grassmann.hpp"
#include "geometry/named_quantities.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"

namespace principia {
namespace numerics {
namespace internal_global_optimization {

using geometry::Position;  // TODO(phl): Use something more general?
using geometry::Vector;
using quantities::Derivative;
using quantities::Length;

template<typename Value, typename Frame>
using Field = std::function<Value(Position<Frame> const&)>;

template<typename Scalar, typename Frame>
using Gradient = Vector<Derivative<Scalar, Length>, Frame>;

template<typename Scalar, typename Frame>
class MultiLevelSingleLinkage {
 public:
  using Box = std::pair<Position<Frame>>;

  MultiLevelSingleLinkage(Box const& box,
                          Field<Scalar, Frame> const& f,
                          Field<Gradient<Scalar, Frame>, Frame> const& grad_f,
                          int64_t values_per_round);
};

}  // namespace internal_global_optimization
}  // namespace numerics
}  // namespace principia

#include "numerics/global_optimization_body.hpp"