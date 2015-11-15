#include "serialization/journal.pb.h"

namespace principia {
namespace ksp_plugin {

class Plugin;

struct DeletePlugin {
  struct In {
    Plugin const* plugin;
  };
  struct Out {};
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
  template<typename Record>
  static void Entry(typename Record::In const& in);
  template<typename Record>
  static void Exit(typename Record::Out const& out);
  template<typename Record>
  static typename Record::Return Return(typename Record::Return const& r3turn);
};

class J0urnal {
 public:
  template<typename Record>
  class Entry {
   public:
    Entry(typename Record::In const& in);
    void Exit(typename Record::Out const& out);
    typename Record::Return Return(typename Record::Return const& r3turn);
  };
};

}  // namespace ksp_plugin
}  // namespace principia
