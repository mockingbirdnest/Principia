#pragma once

#include <string>

namespace principia {
namespace testing_utilities {

// In Debug mode the death message is lost, presumably because Windows tries to
// bring up an abort/retry/ignore dialog.  We don't care too much, the test will
// do the right thing in Release mode.
std::string DeathMessage(std::string const& s);

}  // namespace testing_utilities
}  // namespace principia

#include "testing_utilities/death_message_body.hpp"
