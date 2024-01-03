#pragma once

namespace principia {
namespace numerics {
namespace _transposed_view {
namespace internal {

// TODO(phl): Turn this into a proper view that can be applied to matrices, etc.
template<typename T>
struct TransposedView {
  T const& transpose;
};

template<class T>
TransposedView(T) -> TransposedView<T>;

}  // namespace internal

using internal::TransposedView;

}  // namespace _transposed_view
}  // namespace numerics
}  // namespace principia
