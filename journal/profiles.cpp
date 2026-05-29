#include "journal/profiles.hpp"

#include <cstdint>
#include <cstring>
#include <string>
#include <type_traits>

#include "absl/log/check.h"
#include "absl/log/log.h"

namespace principia {
namespace journal {

#define PRINCIPIA_LAX_POINTER_MAPPING 0
#define PRINCIPIA_PERFORM_RUN_CHECKS 1
#define PRINCIPIA_SET_VERBOSE_LOGGING 1

#if PRINCIPIA_PERFORM_RUN_CHECKS
#define PRINCIPIA_CHECK_EQ(a, b)                                               \
  CHECK((a) == (b)) << "Set PRINCIPIA_PERFORM_RUN_CHECKS to 0 in profile.cpp " \
                    << "to disable this check."
#else
#define PRINCIPIA_CHECK_EQ(a, b)          \
  {                                       \
    [[maybe_unused]] auto const aa = (a); \
    [[maybe_unused]] auto const bb = (b); \
  }
#endif

// The argument `message` below is intentionally not wrapped into parentheses as
// it may not be a complete expression.
#if PRINCIPIA_LAX_POINTER_MAPPING
#define PRINCIPIA_NO_POINTER_MAPPING(message) LOG(ERROR) << message
#else
#define PRINCIPIA_NO_POINTER_MAPPING(message)                                \
  LOG(FATAL) << message                                                      \
             << ".  Set PRINCIPIA_LAX_POINTER_MAPPING to 1 in profiles.cpp " \
                "to disable this check.  "
#endif

namespace {

template<typename T>
void Insert(std::uint64_t const address,
            T* const pointer,
            Player::PointerMap& pointer_map) {
  if (reinterpret_cast<void*>(address) == nullptr) {
    CHECK(pointer == nullptr) << pointer;
    return;
  }
  void* const inserted_pointer = static_cast<void*>(
      const_cast<typename std::remove_cv_t<T>*>(pointer));
  auto const [it, inserted] = pointer_map.emplace(address, inserted_pointer);
  if (!inserted) {
    CHECK_EQ(it->second, inserted_pointer) << address;
  }
}

void Delete(std::uint64_t const address, Player::PointerMap& pointer_map) {
  if (reinterpret_cast<void*>(address) != nullptr) {
    if (auto const it = pointer_map.find(address); it == pointer_map.end()) {
      PRINCIPIA_NO_POINTER_MAPPING(
          "Pointer address not found in Delete: " << address);
    } else {
      pointer_map.erase(it);
    }
  }
}

template<typename T>
  requires(std::is_pointer_v<T>)
T DeserializePointer(std::uint64_t const address,
                     Player::PointerMap const& pointer_map) {
  if (reinterpret_cast<T>(address) == nullptr) {
    return nullptr;
  } else {
    if (auto const it = pointer_map.find(address); it == pointer_map.end()) {
      PRINCIPIA_NO_POINTER_MAPPING(
          "Pointer address not found in DeserializePointer<"
          << typeid(T).name() << ">: " << address);
      return nullptr;
    } else {
      return reinterpret_cast<T>(it->second);
    }
  }
}

// This function uses a `std::string` to store non-UTF-8 data (UTF-16
// specifically) because proto is silly and uses `std::string` for bytes as well
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

// See the comment on `DeserializeUtf16` regarding the usage of `std::string`.
[[maybe_unused]] std::string SerializeUtf16(
    char16_t const* const deserialized) {
  // Note that the string constructed here contains the final char16_t null.
  return std::string(
      reinterpret_cast<char const*>(deserialized),
      sizeof(char16_t) * (std::char_traits<char16_t>::length(deserialized)));
}

}  // namespace

#include "journal/profiles.generated.cc"

#undef PRINCIPIA_LAX_POINTER_MAPPING
#undef PRINCIPIA_PERFORM_RUN_CHECKS
#undef PRINCIPIA_SET_VERBOSE_LOGGING
#undef PRINCIPIA_CHECK_EQ
#undef PRINCIPIA_NO_POINTER_MAPPING

}  // namespace journal
}  // namespace principia
