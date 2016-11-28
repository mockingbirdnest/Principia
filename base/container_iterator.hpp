#pragma once

#include "base/not_null.hpp"

namespace principia {
namespace base {

template<typename Container>
struct ContainerIterator {
  using Iterator = typename Container::iterator;
  ContainerIterator(not_null<Container*> container, Iterator iterator);

  template<typename = std::enable_if_t<
               std::is_same<Iterator,
                            decltype(std::declval<Container>().erase(
                                std::declval<Iterator>()))>::value>>
  ContainerIterator Erase() const;

  not_null<Container*> container;
  typename Container::iterator iterator;
};

}  // namespace base
}  // namespace principia

#include "base/container_iterator_body.hpp"
