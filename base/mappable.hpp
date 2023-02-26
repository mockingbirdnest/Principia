#pragma once

#include "base/not_constructible.hpp"

namespace principia {
namespace base {
namespace _mappable {
namespace internal {

// This struct helps in declaring that a type is "mappable", i.e., that the maps
// declared in principia::geometry can be act on it through the operator().  To
// use it, declare a specialization for the proper |Functor| (a class that must
// have an operator()) and |T| (the class to be made mappable).  The third
// parameter is an enabler and maybe used to restrict specializations with
// SFINAE.
// Any specialization must export two declarations in its public part:
//   o A type named |type|, which is the result of applying |Functor| to |T|.
//   o A static function named |Do|, which applies a functor to a value of |T|
//     and returns a |type|.
// See the comments below for example declarations.
template<typename Functor, typename T, typename = void>
struct Mappable : not_constructible {
  // using type = void;
  // static type Do(Functor const& functor, T const& t);
};

}  // namespace internal

using internal::Mappable;

}  // namespace _mappable
}  // namespace base
}  // namespace principia

namespace principia::base {
using namespace principia::base::_mappable;
}  // namespace principia::base
