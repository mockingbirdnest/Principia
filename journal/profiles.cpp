
#include "journal/profiles.hpp"

#include <fstream>
#include <list>
#include <string>
#include <type_traits>
#include <vector>

#include "base/map_util.hpp"
#include "glog/logging.h"

namespace principia {

using base::FindOrDie;

namespace journal {
namespace {

template<typename T>
void Insert(Player::PointerMap& pointer_map,
            std::uint64_t const address,
            T* const pointer) {
  void* const inserted_pointer = static_cast<void*>(
      const_cast<typename std::remove_cv<T>::type*>(pointer));
  auto const [it, inserted] = pointer_map.emplace(address, inserted_pointer);
  if (!inserted) {
    CHECK_EQ(it->second, inserted_pointer);
  }
}

void Delete(Player::PointerMap& pointer_map,
            std::uint64_t const address) {
  if (reinterpret_cast<void*>(address) != nullptr) {
    auto const it = pointer_map.find(address);
    CHECK(it != pointer_map.end()) << address;
    pointer_map.erase(it);
  }
}

template<typename T,
         typename = typename std::enable_if<std::is_pointer<T>::value>::type>
T DeserializePointer(Player::PointerMap const& pointer_map,
                     std::uint64_t const address) {
  if (reinterpret_cast<T>(address) == nullptr) {
    return nullptr;
  } else {
    return reinterpret_cast<T>(FindOrDie(pointer_map, address));
  }
}

// This function uses a |std::string| to store non-UTF-8 data (UTF-16
// specifically) because proto is silly and uses |std::string| for bytes as well
// as string.
[[maybe_unused]] std::u16string DeserializeUtf16(
    std::string const& serialized) {
  std::u16string result(serialized.size() / sizeof(char16_t), u'\0');
  std::memcpy(result.data(), serialized.c_str(), serialized.size());
  return result;
}

template<typename T>
std::uint64_t SerializePointer(T* t) {
  return reinterpret_cast<std::uint64_t>(t);
}

// See the comment on |DeserializeUtf16| regarding the usage of |std::string|.
[[maybe_unused]] std::string SerializeUtf16(
    char16_t const* const deserialized) {
  // Note that the string constructed here contains the final char16_t null.
  return std::string(
      reinterpret_cast<char const*>(deserialized),
      sizeof(char16_t) *
          (std::char_traits<char16_t const>::length(deserialized)));
}

}  // namespace

#define PRINCIPIA_CHECK_EQ(a, b) CHECK((a) == (b))
#define PRINCIPIA_SET_VERBOSE_LOGGING 1

#include "journal/profiles.generated.cc"

#undef PRINCIPIA_SET_VERBOSE_LOGGING
#undef PRINCIPIA_CHECK_EQ

}  // namespace journal
}  // namespace principia
