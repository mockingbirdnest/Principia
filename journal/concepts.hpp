#pragma once

namespace principia {
namespace journal {
namespace _concepts {
namespace internal {

template<typename P>
concept has_in = requires {
  typename P::In;
};

template<typename P>
concept has_out = requires {
  typename P::Out;
};

template<typename P>
concept has_return = requires {
  typename P::Return;
};

}  // namespace internal

using internal::has_in;
using internal::has_out;
using internal::has_return;

}  // namespace _concepts
}  // namespace journal
}  // namespace principia
