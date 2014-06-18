#pragma once

#include <string>
#include "quantities/quantities.hpp"

namespace principia {
namespace ksp {

std::string Unmanage(System::String^ const string);
System::String^ Manage (std::string const& string);

quantities::Time UniversalTime();

template<typename T>
void Reset(T* pointer, T* const new_pointer) {
  delete(pointer);
  pointer = new_pointer;
}

template<typename T>
void Reset(T* pointer, std::nullptr_t const new_pointer) {
  delete(pointer);
  pointer = new_pointer;
}

}  // namespace ksp
}  // namespace principia
