
#pragma once

#include "base/container_iterator.hpp"

namespace principia {
namespace base {

template<typename Container>
IteratorOn<Container>::IteratorOn(not_null<Container*> container,
                                  Iterator iterator)
    : container_(container),
      iterator_(iterator) {}

template<typename Container>
typename Container::iterator const& IteratorOn<Container>::iterator() const {
  return iterator_;
}

template<typename Container>
template<typename>
IteratorOn<Container> IteratorOn<Container>::Erase() const {
  return IteratorOn(container_, container_->erase(iterator_));
}

}  // namespace base
}  // namespace principia
