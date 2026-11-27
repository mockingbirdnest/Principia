export module principia.base.constant_function;

namespace principia {
namespace base {
namespace _constant_function {
namespace internal {
export {

// A function object that can be called with arbitrary arguments, and always
// returns the same `value`.
template<typename T>
struct ConstantFunction {
  template<typename... Args>
  constexpr T operator()(Args&&...) {
    return value;
  }

  T value;
};

template <typename T>
constexpr ConstantFunction<T> Identically(T value) {
  return {value};
}

}  // export
}  // namespace internal

export using internal::ConstantFunction;
export using internal::Identically;

}  // namespace _constant_function
}  // namespace base
}  // namespace principia
