#pragma once

namespace principia {
namespace geometry {
namespace _concepts {
namespace internal {

template<typename T1, typename T2>
concept hilbert = requires(T1 const& t1, T2 const& t2) {
  InnerProduct(t1, t2);
};

}  // namespace internal

using internal::hilbert;

}  // namespace _concepts
}  // namespace geometry
}  // namespace principia
