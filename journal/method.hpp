
#pragma once

#include <functional>
#include <memory>

#include "base/not_constructible.hpp"
#include "base/not_null.hpp"

namespace principia {
namespace journal {
namespace internal_method {

using base::not_constructible;
using base::not_null;

// The parameter |Profile| is expected to have the following structure:
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

template<typename P, typename = void>
struct has_in : std::false_type, not_constructible {};
template<typename P>
struct has_in<P, std::void_t<typename P::In>> : std::true_type,
                                                not_constructible {};

template<typename P, typename = void>
struct has_out : std::false_type {};
template<typename P>
struct has_out<P, std::void_t<typename P::Out>> : std::true_type,
                                                  not_constructible {};

template<typename P, typename = void>
struct has_return : std::false_type {};
template<typename P>
struct has_return<P, std::void_t<typename P::Return>> : std::true_type,
                                                        not_constructible {};

template<typename Profile>
class Method final {
 public:
  Method();

  // Only declare this constructor if the profile has an |In| type and no |Out|
  // type.
  template<typename P = Profile,
           typename = std::enable_if_t<has_in<P>::value && !has_out<P>::value>>
  explicit Method(typename P::In const& in);

  // Only declare this constructor if the profile has an |Out| type and no |In|
  // type.
  template<typename P = Profile,
           typename = std::enable_if_t<has_out<P>::value && !has_in<P>::value>>
  explicit Method(typename P::Out const& out);

  // Only declare this constructor if the profile has an |In| and an |Out|
  // type.
  template<typename P = Profile,
           typename = std::enable_if_t<has_in<P>::value && has_out<P>::value>>
  Method(typename P::In const& in, typename P::Out const& out);

  ~Method();

  // Only declare this method if the profile has no |Return| type.
  template<typename P = Profile,
           typename = std::enable_if_t<!has_return<P>::value>>
  void Return();

  // Only declare this method if the profile has a |Return| type.
  template<typename P = Profile,
           typename = std::enable_if_t<has_return<P>::value>>
  typename P::Return Return(typename P::Return const& result);

 private:
  std::function<void(not_null<typename Profile::Message*> message)> out_filler_;
  std::function<void(not_null<typename Profile::Message*> message)>
      return_filler_;
  bool returned_ = false;
};

}  // namespace internal_method

using internal_method::Method;

}  // namespace journal
}  // namespace principia

#include "journal/method_body.hpp"
