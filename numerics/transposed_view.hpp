#pragma once

namespace principia {
namespace numerics {

// TODO(phl): Turn this into a proper view that can be applied to matrices, etc.
template<typename T>
struct TransposedView {
  T const& transpose;
};

}  // namespace numerics
}  // namespace principia
