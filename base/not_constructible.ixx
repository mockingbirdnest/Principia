export module principia.base.not_constructible;

namespace principia {
namespace base {
namespace _not_constructible {
namespace internal {

// Type traits should inherit from this, so that they may use inheritance
// without virtual destructors.
struct not_constructible {
  not_constructible() = delete;
  ~not_constructible() = delete;
};

}  // namespace internal

export using internal::not_constructible;

}  // namespace _not_constructible
}  // namespace base
}  // namespace principia
