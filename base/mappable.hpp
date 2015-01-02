#pragma once

#include <functional>
#include <utility>

namespace principia {
namespace base {

//TODO(phl): Comment

template<typename Functor, typename T>
class Mapper {
 public:
  //using type = void;
  //static type Do(Functor const& functor, T const& t);
};

//template<typename T, typename = Mapper<std::function<void(T)>, T>::type>
//class enable_if_mappable2 {
// public:
//  using type = T;
//};
//
//template<typename T>
//class enable_if_mappable {
//public:
//  using type = typename enable_if_mappable2<T>::type;
//};

}  // namespace base
}  // namespace principia