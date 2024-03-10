#pragma once

#include <concepts>

namespace principia {
namespace geometry {
namespace _concepts {
namespace internal {

template<typename T1, typename T2>
concept has_inner_product = requires(T1 const& t1, T2 const& t2) {
  InnerProduct(t1, t2);
};

}  // namespace internal

using internal::has_inner_product;

}  // namespace _concepts
}  // namespace geometry
}  // namespace principia
