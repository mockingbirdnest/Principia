
#pragma once

#include <filesystem>
#include <fstream>
#include <mutex>

#include "base/not_null.hpp"
#include "serialization/journal.pb.h"

namespace principia {
namespace journal {

FORWARD_DECLARE_FROM(method, template<typename Profile> class, Method);

class Recorder final {
 public:
  explicit Recorder(std::filesystem::path const& path);
  ~Recorder();

  // Locking is used to ensure that the pairs of writes don't get intermixed.
  void WriteAtConstruction(serialization::Method const& method);
  void WriteAtDestruction(serialization::Method const& method);

  static void Activate(base::not_null<Recorder*> recorder);
  static void Deactivate();
  static bool IsActivated();

 private:
  void WriteLocked(serialization::Method const& method);

  std::mutex lock_;
  std::ofstream stream_;

  static thread_local Recorder* active_recorder_;

  template<typename>
  friend class Method;
  friend class RecorderTest;
};

}  // namespace journal
}  // namespace principia
