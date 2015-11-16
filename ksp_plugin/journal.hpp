#pragma once

#include <experimental/optional>
#include <functional>
#include <list>
#include <memory>

#include "base/not_null.hpp"
#include "serialization/journal.pb.h"

namespace principia {

using base::not_null;

namespace ksp_plugin {

class Plugin;

struct DeletePlugin {
  struct In {
    Plugin const* plugin;
  };
  struct Out {
    Plugin const** const plugin;
  };

  using Message = serialization::DeletePlugin;
  static void Fill(In const& in, not_null<Message*> const message);
  static void Fill(Out const& out, not_null<Message*> const message);
};

struct NewPlugin {
  struct In {
    double initial_time;
    double planetarium_rotation_in_degrees;
  };
  using Return = Plugin*;

  using Message = serialization::NewPlugin;
  static void Fill(In const& in, not_null<Message*> const message);
  static void Fill(Return const& result, not_null<Message*> const message);
};

class Journal {
 public:
  template<typename Profile>
  class Method {
   public:
    explicit Method(typename Profile::In const& in);

    // Only declare this constructor if the profile has an |Out| type.
    template<typename P = Profile, typename = typename P::Out>
    Method(typename Profile::In const& in, typename P::Out const& out);

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

  template<typename Message>
  static void AppendMessage(not_null<Message*> const message);

 private:
  static std::list<serialization::Method>* journal_;
};

}  // namespace ksp_plugin
}  // namespace principia

#include "ksp_plugin/journal_body.hpp"
