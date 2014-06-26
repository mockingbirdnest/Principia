#include "test_plugin/test_plugin.hpp"

#include <string>

namespace principia {
namespace test_plugin {

int Say33() {
  return 33;
}

char const* SayHello() {
  return "Hello from native C++!";
}

}  // namespace test_plugin
}  // namespace princpia
