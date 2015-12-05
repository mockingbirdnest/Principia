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

template<typename Profile>
class Method {
 public:
  Method();

  // Only declare this constructor if the profile has an |In| type.
  template<typename P = Profile, typename = typename P::In>
  explicit Method(typename P::In const& in);

  // Only declare this constructor if the profile has an |In| and an |Out|
  // type.
  template<typename P = Profile,
            typename = typename P::In, typename = typename P::Out>
  Method(typename P::In const& in, typename P::Out const& out);

  ~Method();

  void Return();

  // Only declare this method if the profile has a |Return| type.
  template<typename P = Profile, typename = typename P::Return>
  typename P::Return Return(typename P::Return const& result);

 private:
  std::unique_ptr<typename Profile::Message> message_;
  std::function<void()> out_filler_;
  bool returned_ = false;
};

}  // namespace journal
}  // namespace principia

#include "journal/method_body.hpp"
