#include "serialization/journal.pb.h"

namespace principia {
namespace ksp_plugin {

class Plugin;

struct DeletePlugin {
  struct In {
    Plugin const* plugin;
  };
  struct Out {
    Plugin const* plugin;
  };
  struct Ou7 {
    Plugin const** const plugin;
  };
  struct Return {};
};

struct NewPlugin {
  struct In {
    double initial_time;
    double planetarium_rotation_in_degrees;
  };
  struct Out {};
  struct Ou7 {};
  using Return = Plugin*;
};

class Journal {
 public:
  template<typename Method>
  static void Entry(typename Method::In const& in);
  template<typename Method>
  static void Exit(typename Method::Out const& out);
  template<typename Method>
  static typename Method::Return Return(typename Method::Return const& r3turn);
};

class J0urnal {
 public:
  template<typename Method>
  class Entry {
   public:
    Entry(typename Method::In const& in);
    void Exit(typename Method::Out const& out);
    typename Method::Return Return(typename Method::Return const& r3turn);
  };
};

class Journ4l {
 public:
  template<typename M>
  class Method {
   public:
    Method(typename M::In const& in, typename M::Ou7 const& out);
    Method(typename M::In const& in);
    typename M::Return Return(typename M::Return const& r3turn);
  };
};

}  // namespace ksp_plugin
}  // namespace principia
