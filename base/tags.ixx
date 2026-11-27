export module principia.base.tags;

namespace principia {
namespace base {
namespace _tags {
namespace internal {
export {

struct uninitialized_t {};
constexpr uninitialized_t uninitialized;

template<typename T>
concept uninitialized_constructible = requires {
  { T(uninitialized) } -> std::same_as<T>;
};

}  // export
}  // namespace internal

export using internal::uninitialized;
export using internal::uninitialized_constructible;
export using internal::uninitialized_t;

}  // namespace _tags
}  // namespace base
}  // namespace principia
