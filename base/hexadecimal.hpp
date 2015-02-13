#pragma once

#include <string>
#include <type_traits>

#include "base/not_null.hpp"

namespace principia {
namespace base {

// The encoding can be done in-place.  The result is upper-case.
template<typename Container>
std::enable_if_t<
    std::is_convertible<typename Container::value_type, uint8_t>::value,
    void>
HexadecimalEncode(Container const& input, not_null<Container*> output);

// Invalid digit pairs are decoded to NUL.  The decoding can be done in-place.
// Returns true if, and only if, all digit pairs were valid and |input.size()|
// was even. If |input.size()| is odd, the last character of the input is
// ignored.  Ignores case.
template<typename Container>
std::enable_if_t<
    std::is_convertible<typename Container::value_type, uint8_t>::value,
    bool>
HexadecimalDecode(Container const& input, not_null<Container*> output);

}  // namespace base
}  // namespace principia

#include "base/hexadecimal_body.hpp"
