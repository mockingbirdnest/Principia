#pragma once

#include <fstream>
#include <functional>
#include <list>
#include <memory>

#include "base/not_null.hpp"
#include "ksp_plugin/interface.hpp"
#include "serialization/journal.pb.h"

namespace principia {

using base::not_null;

namespace ksp_plugin {

struct InitGoogleLogging {
  using Message = serialization::InitGoogleLogging;
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

struct DirectlyInsertCelestial {
  struct In {
    Plugin* plugin;
    int celestial_index;
    int const* parent_index;
    char const* gravitational_parameter;
    char const* axis_right_ascension;
    char const* axis_declination;
    char const* j2;
    char const* reference_radius;
    char const* x;
    char const* y;
    char const* z;
    char const* vx;
    char const* vy;
    char const* vz;
  };

  using Message = serialization::DirectlyInsertCelestial;
  static void Fill(In const& in, not_null<Message*> const message);
};

struct InsertCelestial {
  struct In {
    Plugin* plugin;
    int celestial_index;
    double gravitational_parameter;
    int parent_index;
    QP from_parent;
  };

  using Message = serialization::InsertCelestial;
  static void Fill(In const& in, not_null<Message*> const message);
};

class Journal {
 public:
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

  Journal(std::string const& filename);
  ~Journal();

  void Write(serialization::Method const& method);

  static void Activate(base::not_null<Journal*> const journal);
  static void Deactivate();

 private:
  std::ofstream stream_;

  static Journal* active_;

  template<typename>
  friend class Method;
  friend class JournalTest;
};

}  // namespace ksp_plugin
}  // namespace principia

#include "ksp_plugin/journal_body.hpp"
