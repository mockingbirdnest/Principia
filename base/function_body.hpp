#pragma once

#include <utility>

#include "base/function.hpp"

namespace principia {
namespace base {
namespace internal_function {

template<typename F, typename Result, typename... Args>
ConcreteFunctor<F, Result, Args...>::ConcreteFunctor(F Functor)
    : functor_(std::move(Functor)) {}

template<typename F, typename Result, typename... Args>
Result ConcreteFunctor<F, Result, Args...>::call(Args&&... args) {
  return functor_(std::forward<Args>(args)...);
}

template<typename Result, typename... Args>
template<typename F>
function<Result(Args...)>::function(F Functor)
    : functor_(std::make_unique<ConcreteFunctor<F, Result, Args...>>(
          std::move(Functor))) {}

template<typename Result, typename... Args>
Result function<Result(Args...)>::operator()(Args&&... args) {
  return functor_->call(std::forward<Args>(args)...);
}

}  // namespace internal_functio
}  // namespace base
}  // namespace principia
