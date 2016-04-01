
#pragma once

#include <functional>
#include <memory>

namespace principia {
namespace journal {

// The parameter |Profile| is expected to have the following structure:
//
//  struct SomeProfile {
//    struct In {  // Only present if the profile has 'in' parameters.
//      ...        // 'in' parameters copied verbatim from the actual profile.
//    };
//    struct Out {  // Only present if the profile has 'out' parameters.
//      ...         // 'out' parameters copied verbatim from the actual profile.
//    };
//    using Return = ...;  // Only present if the profile does not return void.
//
//    using Message = serialization::SerializePlugin;
//
//    // The following functions must be omitted if In/Out/Return is omitted.
//    static void Fill(In const& in, not_null<Message*> const message);
//    static void Fill(Out const& out, not_null<Message*> const message);
//    static void Fill(Return const& result, not_null<Message*> const message);
//
//    static void Run(Message const& message,
//                    not_null<Player::PointerMap*> const pointer_map);
//  };

namespace internal {

template<typename T>
using void_if_exists = void;

template<typename P, typename = void>
struct has_in : std::false_type {};
template<typename P>
struct has_in<P, void_if_exists<typename P::In>> : std::true_type {};

template<typename P, typename = void>
struct has_out : std::false_type {};
template<typename P>
struct has_out<P, void_if_exists<typename P::Out>> : std::true_type {};

template<typename P, typename = void>
struct has_return : std::false_type {};
template<typename P>
struct has_return<P, void_if_exists<typename P::Return>> : std::true_type {};

}  // namespace internal

template<typename Profile>
class Method {
 public:
  Method();

  // Only declare this constructor if the profile has an |In| type and no |Out|
  // type.
  template<typename P = Profile,
           typename = std::enable_if_t<internal::has_in<P>::value &&
                                       !internal::has_out<P>::value>>
  explicit Method(typename P::In const& in);

  // Only declare this constructor if the profile has an |Out| type and no |In|
  // type.
  template<typename P = Profile,
           typename = std::enable_if_t<internal::has_out<P>::value &&
                                       !internal::has_in<P>::value>>
  explicit Method(typename P::Out const& out);

  // Only declare this constructor if the profile has an |In| and an |Out|
  // type.
  template<typename P = Profile,
           typename = std::enable_if_t<internal::has_in<P>::value &&
                                       internal::has_out<P>::value>>
  Method(typename P::In const& in, typename P::Out const& out);

  ~Method();

  // Only declare this method if the profile has no |Return| type.
  template<typename P = Profile,
           typename = std::enable_if_t<!internal::has_return<P>::value>>
  void Return();

  // Only declare this method if the profile has a |Return| type.
  template<typename P = Profile,
           typename = std::enable_if_t<internal::has_return<P>::value>>
  typename P::Return Return(typename P::Return const& result);

 private:
  void LogMethodIfDebug();

  std::unique_ptr<typename Profile::Message> message_;
  std::function<void()> out_filler_;
  bool returned_ = false;
};

}  // namespace journal
}  // namespace principia

#include "journal/method_body.hpp"
