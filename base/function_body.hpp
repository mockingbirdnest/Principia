#pragma once

#include <utility>

#include "base/function.hpp"

namespace principia {
namespace base {
namespace _function {
namespace internal {

template<typename F, typename Result, typename... Args>
ConcreteFunctor<F, Result, Args...>::ConcreteFunctor(F functor)
    : functor_(std::move(functor)) {}

template<typename F, typename Result, typename... Args>
Result ConcreteFunctor<F, Result, Args...>::Call(Args&&... args) {
  return functor_(std::forward<Args>(args)...);
}

template<typename Result, typename... Args>
function<Result(Args...)>::function() {}

template<typename Result, typename... Args>
template<typename F>
function<Result(Args...)>::function(F functor)
    : functor_(std::make_unique<ConcreteFunctor<F, Result, Args...>>(
          std::move(functor))) {}

template<typename Result, typename... Args>
Result function<Result(Args...)>::operator()(Args&&... args) const {
  return functor_->Call(std::forward<Args>(args)...);
}

}  // namespace internal
}  // namespace _function
}  // namespace base
}  // namespace principia
