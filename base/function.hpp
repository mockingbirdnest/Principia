#pragma once

#include <memory>

namespace principia {
namespace base {
namespace internal_function {

template<typename Result, typename... Args>
class Functor {
 public:
  virtual ~Functor() = default;
  virtual Result call(Args&&...) = 0;
};

template<typename F, typename Result, typename... Args>
class ConcreteFunctor final : public Functor<Result, Args...> {
 public:
  ConcreteFunctor(F Functor);
  Result call(Args&&... args) final;
 private:
  F functor_;
};

template<typename Signature>
struct function;

template<typename Result, typename... Args>
class function<Result(Args...)> {
 public:
  template<typename F>
  function(F functor);

  function(function&&) = default;
  function& operator=(function&&) = default;

  Result operator()(Args&&... args);

 private:
  std::unique_ptr<Functor<Result, Args...>> functor_;
};

}  // namespace internal_function

using internal_function::function;

}  // namespace base
}  // namespace principia

#include "base/function_body.hpp"
