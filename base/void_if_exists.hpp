#pragma once

namespace principia {
namespace base {

// This alias allows SFINAE on the existence of a type.  It can be used together
// with decltype to perform SFINAE on the legality of expressions.
template<typename T>
using void_if_exists = void;

}  // namespace base
}  // namespace principia
