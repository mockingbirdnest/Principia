#pragma once

#include "base/not_null.hpp"

namespace principia {
namespace base {

template<typename Container>
class IteratorOn final {
 public:
  using Iterator = typename Container::iterator;

  IteratorOn(not_null<Container*> container, Iterator iterator);

  template<typename = std::enable_if_t<
               std::is_same<Iterator,
                            decltype(std::declval<Container>().erase(
                                std::declval<Iterator>()))>::value>>
  IteratorOn Erase() const;

  Iterator const& iterator() const;

 private:
  not_null<Container*> container_;
  Iterator iterator_;
};

}  // namespace base
}  // namespace principia

#include "base/container_iterator_body.hpp"
