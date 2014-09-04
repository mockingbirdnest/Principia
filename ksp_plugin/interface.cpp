#include "ksp_plugin/interface.hpp"

#include <string>

namespace principia {
namespace ksp_plugin {

int Say33() {
  return 33;
}

char const* SayHello() {
  return "Hello from native C++!";
}

}  // namespace ksp_plugin
}  // namespace principia
