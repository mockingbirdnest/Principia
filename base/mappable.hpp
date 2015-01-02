#pragma once

#include <functional>
#include <utility>

namespace principia {
namespace base {

//TODO(phl): Comment

template<typename Functor, typename T, typename = void>
class Mappable {
 public:
  //using type = void;
  //static type Do(Functor const& functor, T const& t);
};

}  // namespace base
}  // namespace principia