
#include "ksp/utilities.hpp"

#include <string>

#include <msclr/marshal.h>

using principia::quantities::SIUnit;
using principia::quantities::Time;

namespace principia {
namespace ksp {

std::string Unmanage(System::String^ const string) {
  msclr::interop::marshal_context^ context =
      gcnew msclr::interop::marshal_context();
  std::string result = context->marshal_as<const char*>(string);
  delete context;
  return result;
}

System::String^ Manage (std::string const& string) {
  return gcnew System::String(string.c_str());
}

Time UniversalTime() {
  return Planetarium::GetUniversalTime() * SIUnit<Time>();
}

}  // namespace ksp
}  // namespace principia
