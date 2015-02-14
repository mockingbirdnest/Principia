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

// Invalid digits are read as 0.  If |input.size()| is odd, the last
// character of the input is ignored.  Ignores case.
// The decoding can be done in-place, but if the containers overlap,
// |output.data()| shall not exceed |&input.data()[1]|.
template<typename Container>
std::enable_if_t<
    std::is_convertible<typename Container::value_type, uint8_t>::value>
HexadecimalDecode(Container const& input, not_null<Container*> output);

}  // namespace base
}  // namespace principia

#include "base/hexadecimal_body.hpp"
