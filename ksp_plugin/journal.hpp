#pragma once

#include <experimental/optional>
#include <functional>
#include <memory>

#include "base/not_null.hpp"
#include "serialization/journal.pb.h"

namespace principia {

using base::not_null;

namespace ksp_plugin {

class Plugin;

struct DeletePlugin {
  using Message = serialization::DeletePlugin;
  struct In {
    Plugin const* plugin;
  };
  struct Out {
    Plugin const** const plugin;
  };
  struct Return {};

  static void Fill(In const& in, not_null<Message*> const message);
  static void Fill(Out const& out, not_null<Message*> const message);
  static void Fill(Return const& r3turn, not_null<Message*> const message);
};

struct NewPlugin {
  using Message = serialization::NewPlugin;
  struct In {
    double initial_time;
    double planetarium_rotation_in_degrees;
  };
  struct Out {};
  using Return = Plugin*;

  static void Fill(In const& in, not_null<Message*> const message);
  static void Fill(Out const& out, not_null<Message*> const message);
  static void Fill(Return const& r3turn, not_null<Message*> const message);
};

class Journal {
 public:
  template<typename Profile>
  class Method {
   public:
    Method(typename Profile::In const& in, typename Profile::Out const& out);
    explicit Method(typename Profile::In const& in);
    ~Method();

    typename Profile::Return Return(typename Profile::Return const& r3turn);
    void Return();

   private:
    std::unique_ptr<typename Profile::Message> const message_;
    std::experimental::optional<typename Profile::Out> out_;
    bool returned_ = false;
  };
};

}  // namespace ksp_plugin
}  // namespace principia

#include "ksp_plugin/journal_body.hpp"
