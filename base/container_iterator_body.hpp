
#pragma once

#include "base/container_iterator.hpp"

namespace principia {
namespace base {

template<typename Container>
ContainerIterator<Container>::ContainerIterator(not_null<Container*> container,
                                                Iterator iterator)
    : container(container),
      iterator(iterator) {}

template<typename Container>
template<typename>
ContainerIterator<Container> ContainerIterator<Container>::Erase() const {
  return ContainerIterator(container, container->erase(iterator));
}

}  // namespace base
}  // namespace principia
