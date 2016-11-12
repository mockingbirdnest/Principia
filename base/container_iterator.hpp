#pragma once

#include "base/not_null.hpp"

namespace principia {
namespace base {

template<typename Container>
struct ContainerIterator {
  ContainerIterator(not_null<Container*> container,
                    typename Container::iterator iterator)
      : container(container),
        iterator(iterator) {}

  not_null<Container*> container;
  typename Container::iterator iterator;
};

}  // namespace base
}  // namespace principia
