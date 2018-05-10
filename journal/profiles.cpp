
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
  auto inserted = pointer_map.emplace(address, inserted_pointer);
  if (!inserted.second) {
    CHECK_EQ(inserted.first->second, inserted_pointer);
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

char16_t const* DeserializeUtf16(std::string const& serialized) {
  // The result is char16_t-null-terminated because of the way the serialized
  // string was built in SerializeUtf16.
  return reinterpret_cast<char16_t const*>(serialized.c_str());
}

template<typename T>
std::uint64_t SerializePointer(T* t) {
  return reinterpret_cast<std::uint64_t>(t);
}

std::string SerializeUtf16(char16_t const* const deserialized) {
  // Note that the string constructed here contains the final char16_t null.
  return std::string(
      reinterpret_cast<char const*>(deserialized),
      sizeof(char16_t) *
          (std::char_traits<char16_t const>::length(deserialized) + 1));
}

}  // namespace

#define PRINCIPIA_CHECK_EQ(a, b) CHECK((a) == (b))
#define PRINCIPIA_SET_VERBOSE_LOGGING 1

#include "journal/profiles.generated.cc"

#undef PRINCIPIA_SET_VERBOSE_LOGGING
#undef PRINCIPIA_CHECK_EQ

}  // namespace journal
}  // namespace principia
