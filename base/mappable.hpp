#pragma once

#include <utility>

namespace principia {
namespace base {

//TODO(phl): Comment

template<typename Functor, typename T>
class Mapper {
 public:
  using type = T;
  static type Do(Functor const& functor, T const& t);
};

template<typename T, typename = decltype(Mapper<int, T>::Do(0, std::declval<T>()))>
class enable_if_mappable2 {
 public:
  using type = T;
};

template<typename T>
class enable_if_mappable {
public:
  using type = typename enable_if_mappable2<T>::type;
};

}  // namespace base
}  // namespace principia