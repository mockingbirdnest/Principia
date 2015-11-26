#pragma once

#include <experimental/filesystem>
#include <fstream>
#include <functional>
#include <map>
#include <memory>
#include <string>

#include "base/not_null.hpp"
#include "ksp_plugin/interface.hpp"
#include "serialization/journal.pb.h"

namespace principia {

using base::not_null;

namespace ksp_plugin {

using PointerMap = std::map<std::uint64_t, void*>;

struct InitGoogleLogging {
  using Message = serialization::InitGoogleLogging;
};

struct SetBufferedLogging {
  struct In {
    int const max_severity;
  };

  using Message = serialization::SetBufferedLogging;
  static void Fill(In const& in, not_null<Message*> const message);
  static void Run(Message const& message,
                  not_null<PointerMap*> const pointer_map);
};

struct GetBufferedLogging {
  using Return = int;

  using Message = serialization::GetBufferedLogging;
  static void Fill(Return const& result, not_null<Message*> const message);
  static void Run(Message const& message,
                  not_null<PointerMap*> const pointer_map);
};

struct SetBufferDuration {
  struct In {
    int const seconds;
  };

  using Message = serialization::SetBufferDuration;
  static void Fill(In const& in, not_null<Message*> const message);
  static void Run(Message const& message,
                  not_null<PointerMap*> const pointer_map);
};

struct GetBufferDuration {
  using Return = int;

  using Message = serialization::GetBufferDuration;
  static void Fill(Return const& result, not_null<Message*> const message);
  static void Run(Message const& message,
                  not_null<PointerMap*> const pointer_map);
};

struct SetSuppressedLogging {
  struct In {
    int const min_severity;
  };

  using Message = serialization::SetSuppressedLogging;
  static void Fill(In const& in, not_null<Message*> const message);
  static void Run(Message const& message,
                  not_null<PointerMap*> const pointer_map);
};

struct GetSuppressedLogging {
  using Return = int;

  using Message = serialization::GetSuppressedLogging;
  static void Fill(Return const& result, not_null<Message*> const message);
  static void Run(Message const& message,
                  not_null<PointerMap*> const pointer_map);
};

struct SetVerboseLogging {
  struct In {
    int const level;
  };

  using Message = serialization::SetVerboseLogging;
  static void Fill(In const& in, not_null<Message*> const message);
  static void Run(Message const& message,
                  not_null<PointerMap*> const pointer_map);
};

struct GetVerboseLogging {
  using Return = int;

  using Message = serialization::GetVerboseLogging;
  static void Fill(Return const& result, not_null<Message*> const message);
  static void Run(Message const& message,
                  not_null<PointerMap*> const pointer_map);
};

struct SetStderrLogging {
  struct In {
    int const min_severity;
  };

  using Message = serialization::SetStderrLogging;
  static void Fill(In const& in, not_null<Message*> const message);
  static void Run(Message const& message,
                  not_null<PointerMap*> const pointer_map);
};

struct GetStderrLogging {
  using Return = int;

  using Message = serialization::GetStderrLogging;
  static void Fill(Return const& result, not_null<Message*> const message);
  static void Run(Message const& message,
                  not_null<PointerMap*> const pointer_map);
};

struct LogInfo {
  struct In {
    char const* const message;
  };

  using Message = serialization::LogInfo;
  static void Fill(In const& in, not_null<Message*> const message);
  static void Run(Message const& message,
                  not_null<PointerMap*> const pointer_map);
};

struct LogWarning {
  struct In {
    char const* const message;
  };

  using Message = serialization::LogWarning;
  static void Fill(In const& in, not_null<Message*> const message);
  static void Run(Message const& message,
                  not_null<PointerMap*> const pointer_map);
};

struct LogError {
  struct In {
    char const* const message;
  };

  using Message = serialization::LogError;
  static void Fill(In const& in, not_null<Message*> const message);
  static void Run(Message const& message,
                  not_null<PointerMap*> const pointer_map);
};

struct LogFatal {
  struct In {
    char const* const message;
  };

  using Message = serialization::LogFatal;
  static void Fill(In const& in, not_null<Message*> const message);
  static void Run(Message const& message,
                  not_null<PointerMap*> const pointer_map);
};

struct NewPlugin {
  struct In {
    double const initial_time;
    double const planetarium_rotation_in_degrees;
  };
  using Return = Plugin*;

  using Message = serialization::NewPlugin;
  static void Fill(In const& in, not_null<Message*> const message);
  static void Fill(Return const& result, not_null<Message*> const message);
  static void Run(Message const& message,
                  not_null<PointerMap*> const pointer_map);
};

struct DeletePlugin {
  struct In {
    Plugin const* const plugin;
  };
  struct Out {
    Plugin const** const plugin;
  };

  using Message = serialization::DeletePlugin;
  static void Fill(In const& in, not_null<Message*> const message);
  static void Fill(Out const& out, not_null<Message*> const message);
  static void Run(Message const& message,
                  not_null<PointerMap*> const pointer_map);
};

struct DirectlyInsertCelestial {
  struct In {
    Plugin* const plugin;
    int const celestial_index;
    int const* const parent_index;
    char const* const gravitational_parameter;
    char const* const axis_right_ascension;
    char const* const axis_declination;
    char const* const j2;
    char const* const reference_radius;
    char const* const x;
    char const* const y;
    char const* const z;
    char const* const vx;
    char const* const vy;
    char const* const vz;
  };

  using Message = serialization::DirectlyInsertCelestial;
  static void Fill(In const& in, not_null<Message*> const message);
  static void Run(Message const& message,
                  not_null<PointerMap*> const pointer_map);
};

struct InsertCelestial {
  struct In {
    Plugin* const plugin;
    int const celestial_index;
    double const gravitational_parameter;
    int const parent_index;
    QP const from_parent;
  };

  using Message = serialization::InsertCelestial;
  static void Fill(In const& in, not_null<Message*> const message);
  static void Run(Message const& message,
                  not_null<PointerMap*> const pointer_map);
};

struct InsertSun {
  struct In {
    Plugin* const plugin;
    int const celestial_index;
    double const gravitational_parameter;
  };

  using Message = serialization::InsertSun;
  static void Fill(In const& in, not_null<Message*> const message);
  static void Run(Message const& message,
                  not_null<PointerMap*> const pointer_map);
};

struct UpdateCelestialHierarchy {
  struct In {
    Plugin const* const plugin;
    int const celestial_index;
    int const parent_index;
  };

  using Message = serialization::UpdateCelestialHierarchy;
  static void Fill(In const& in, not_null<Message*> const message);
  static void Run(Message const& message,
                  not_null<PointerMap*> const pointer_map);
};

struct EndInitialization {
  struct In {
    Plugin const* const plugin;
  };

  using Message = serialization::EndInitialization;
  static void Fill(In const& in, not_null<Message*> const message);
  static void Run(Message const& message,
                  not_null<PointerMap*> const pointer_map);
};

struct InsertOrKeepVessel {
  struct In {
    Plugin* const plugin;
    char const* const vessel_guid;
    int const parent_index;
  };
  using Return = bool;

  using Message = serialization::InsertOrKeepVessel;
  static void Fill(In const& in, not_null<Message*> const message);
  static void Fill(Return const& result, not_null<Message*> const message);
  static void Run(Message const& message,
                  not_null<PointerMap*> const pointer_map);
};

struct SetVesselStateOffset {
  struct In {
    Plugin* const plugin;
    char const* const vessel_guid;
    QP const from_parent;
  };

  using Message = serialization::SetVesselStateOffset;
  static void Fill(In const& in, not_null<Message*> const message);
  static void Run(Message const& message,
                  not_null<PointerMap*> const pointer_map);
};

struct AdvanceTime {
  struct In {
    Plugin* const plugin;
    double const t;
    double const planetarium_rotation;
  };

  using Message = serialization::AdvanceTime;
  static void Fill(In const& in, not_null<Message*> const message);
  static void Run(Message const& message,
                  not_null<PointerMap*> const pointer_map);
};

struct ForgetAllHistoriesBefore {
  struct In {
    Plugin* const plugin;
    double const t;
  };

  using Message = serialization::ForgetAllHistoriesBefore;
  static void Fill(In const& in, not_null<Message*> const message);
  static void Run(Message const& message,
                  not_null<PointerMap*> const pointer_map);
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

  explicit Journal(std::experimental::filesystem::path const& path);
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

class Player {
 public:
  explicit Player(std::experimental::filesystem::path const& path);

  // Replays the next message in the journal.  Returns false at end of journal.
  bool Play();

 private:
  // Reads one message from the stream.  Returns a |nullptr| at end of stream.
  std::unique_ptr<serialization::Method> Read();

  template<typename Profile>
  bool RunIfAppropriate(serialization::Method const& method);

  PointerMap pointer_map_;
  std::ifstream stream_;

  friend class JournalTest;
};

}  // namespace ksp_plugin
}  // namespace principia

#include "ksp_plugin/journal_body.hpp"
