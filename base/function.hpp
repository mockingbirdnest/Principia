#pragma once

#include <memory>

namespace principia {
namespace base {
namespace _function {
namespace internal {

template<typename Result, typename... Args>
class Functor {
 public:
  virtual ~Functor() = default;
  virtual Result Call(Args&&...) = 0;
};

template<typename F, typename Result, typename... Args>
class ConcreteFunctor final : public Functor<Result, Args...> {
 public:
  explicit ConcreteFunctor(F functor);
  Result Call(Args&&... args) final;
 private:
  F functor_;
};

template<typename Signature>
class function;

template<typename Result, typename... Args>
class function<Result(Args...)> {
 public:
  function();

  template<typename F>
  function(F functor);

  function(function&&) = default;  // NOLINT(whitespace/operators)
  function& operator=(function&&) = default;

  Result operator()(Args&&... args) const;

 private:
  std::unique_ptr<Functor<Result, Args...>> functor_;
};

}  // namespace internal

using internal::function;

}  // namespace _function
}  // namespace base
}  // namespace principia

namespace principia::base {
using namespace principia::base::_function;
}  // namespace principia::base

#include "base/function_body.hpp"
