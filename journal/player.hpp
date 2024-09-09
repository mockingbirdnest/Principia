#pragma once

#include <filesystem>
#include <fstream>
#include <map>
#include <memory>

#include "serialization/journal.pb.h"

namespace principia {
namespace journal {

class PlayerTest;
class RecorderTest;

namespace _player {
namespace internal {

class Player final {
 public:
  using PointerMap = std::map<std::uint64_t, void*>;

  explicit Player(std::filesystem::path const& path);

  // Replays the next message in the journal.  Returns false at end of journal.
  // `index` is the 0-based index of the message in the journal.
  bool Play(int index);

  // Same as `Play`, but does not execute the messages, only parse them.
  bool Scan(int index);

  // Return the last replayed messages.
  serialization::Method const& last_method_in() const;
  serialization::Method const& last_method_out_return() const;

 private:
  // Reads one message from the stream.  Returns a `nullptr` at end of stream.
  std::unique_ptr<serialization::Method> Read();

  // Implementation of `Play` and `Scan`.
  bool Process(std::unique_ptr<serialization::Method> method_in,
               int const index, bool const play);

  template<typename Profile>
  bool RunIfAppropriate(serialization::Method const& method_in,
                        serialization::Method const& method_out_return);

  PointerMap pointer_map_;
  std::ifstream stream_;

  std::unique_ptr<serialization::Method> last_method_in_;
  std::unique_ptr<serialization::Method> last_method_out_return_;

  friend class journal::PlayerTest;
  friend class journal::RecorderTest;
};

}  // namespace internal

using internal::Player;

}  // namespace _player
}  // namespace journal
}  // namespace principia

#include "journal/player_body.hpp"
