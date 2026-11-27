export module principia.base.concepts;

namespace principia {
namespace base {
namespace _concepts {
namespace internal {
export {

// True if and only if T has a (possibly templated) static member function named
// ReadFromMessage.
template<typename T>
concept serializable = requires {
  &T::ReadFromMessage;
} || requires {
  &T::template ReadFromMessage<>;
};

}  // export
}  // namespace internal

export using internal::serializable;

}  // namespace _concepts
}  // namespace base
}  // namespace principia
