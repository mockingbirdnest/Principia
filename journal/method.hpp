#pragma once

#include <functional>
#include <memory>
#include <type_traits>

#include "base/not_constructible.hpp"
#include "base/not_null.hpp"
#include "journal/concepts.hpp"

namespace principia {
namespace journal {
namespace _method {
namespace internal {

using namespace principia::base::_not_constructible;
using namespace principia::base::_not_null;
using namespace principia::journal::_concepts;

// The parameter `Profile` is expected to have the following structure:
//
//  struct SomeProfile : not_constructible {
//    struct In final {  // Only present if the profile has 'in' parameters.
//      ...  // 'in' parameters copied verbatim from the actual profile.
//    };
//    struct Out final {  // Only present if the profile has 'out' parameters.
//      ...  // 'out' parameters copied verbatim from the actual profile.
//    };
//    using Return = ...;  // Only present if the profile does not return void.
//
//    using Message = serialization::SerializePlugin;
//
//    // The following functions must be omitted if In/Out/Return is omitted.
//    static void Fill(In const& in, not_null<Message*> message);
//    static void Fill(Out const& out, not_null<Message*> message);
//    static void Fill(Return const& result, not_null<Message*> message);
//
//    static void Run(Message const& message,
//                    not_null<Player::PointerMap*> pointer_map);
//  };

template<typename Profile>
class Method final {
 public:
  Method();

  // Only declare this constructor if the profile has an `In` type and no `Out`
  // type.
  template<typename P = Profile>
  explicit Method(typename P::In const& in)
    requires has_in<P> && (!has_out<P>);

  // Only declare this constructor if the profile has an `Out` type and no `In`
  // type.
  template<typename P = Profile>
  explicit Method(typename P::Out const& out)
    requires has_out<P> && (!has_in<P>);

  // Only declare this constructor if the profile has an `In` and an `Out`
  // type.
  template<typename P = Profile>
  Method(typename P::In const& in, typename P::Out const& out)
    requires has_in<P> && has_out<P>;

  ~Method();

  // Only declare this method if the profile has no `Return` type.
  template<typename P = Profile>
  void Return() requires (!has_return<P>);  // NOLINT

  // Only declare this method if the profile has a `Return` type.
  template<typename P = Profile>
  typename P::Return Return(typename P::Return const& result)
    requires has_return<P>;

 private:
  std::function<void(not_null<typename Profile::Message*> message)> out_filler_;
  std::function<void(not_null<typename Profile::Message*> message)>
      return_filler_;
  bool returned_ = false;
};

}  // namespace internal
}  // namespace _method

// To preserve a reasonable style in the interface, we export this class
// directly in `journal`.
using _method::internal::Method;

}  // namespace journal
}  // namespace principia

#include "journal/method_body.hpp"
