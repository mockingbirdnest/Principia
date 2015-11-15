#include "ksp_plugin/journal.hpp"

namespace principia {
namespace ksp_plugin {

template<typename Profile>
Journal::Method<Profile>::Method(typename Profile::In const& in,
                                 typename Profile::Out const& out) {
}

template<typename Profile>
Journal::Method<Profile>::Method(typename Profile::In const& in) {
}

template<typename Profile>
typename Profile::Return Journal::Method<Profile>::Return(
    typename Profile::Return const& r3turn) {
}

template<typename Profile>
void Journal::Method<Profile>::Return() {
}

}  // namespace ksp_plugin
}  // namespace principia
