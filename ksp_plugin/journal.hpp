#include "serialization/journal.pb.h"

namespace principia {
namespace ksp_plugin {

class Plugin;

struct DeletePlugin {
  struct In {
    Plugin const* plugin;
  };
  struct Out {
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
  using Return = Plugin*;
};

class Journal {
 public:
  template<typename Profile>
  class Method {
   public:
    Method(typename Profile::In const& in, typename Profile::Out const& out);
    Method(typename Profile::In const& in);
    typename Profile::Return Return(typename Profile::Return const& r3turn);
    void Return();
  };
};

}  // namespace ksp_plugin
}  // namespace principia

#include "ksp_plugin/journal_body.hpp"
