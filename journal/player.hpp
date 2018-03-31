
#pragma once

#include <filesystem>
#include <fstream>
#include <map>
#include <memory>

#include "serialization/journal.pb.h"

namespace principia {
namespace journal {

class Player final {
 public:
  using PointerMap = std::map<std::uint64_t, void*>;

  explicit Player(std::filesystem::path const& path);

  // Replays the next message in the journal.  Returns false at end of journal.
  bool Play();

  // Return the last replayed messages.
  serialization::Method const& last_method_in() const;
  serialization::Method const& last_method_out_return() const;

 private:
  // Reads one message from the stream.  Returns a |nullptr| at end of stream.
  std::unique_ptr<serialization::Method> Read();

  template<typename Profile>
  bool RunIfAppropriate(serialization::Method const& method_in,
                        serialization::Method const& method_out_return);

  PointerMap pointer_map_;
  std::ifstream stream_;

  std::unique_ptr<serialization::Method> last_method_in_;
  std::unique_ptr<serialization::Method> last_method_out_return_;

  friend class PlayerTest;
  friend class RecorderTest;
};

}  // namespace journal
}  // namespace principia

#include "journal/player_body.hpp"
