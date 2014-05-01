#include "stdafx.hpp"

#include "TestUtilities.hpp"

namespace Principia {
namespace TestUtilities {
// The Microsoft equivalent only takes a wchar_t*.
void WriteLog(std::wstring const& message) {
  Microsoft::VisualStudio::CppUnitTestFramework::Logger::WriteMessage(
      message.c_str());
}
void NewLine() {
  Microsoft::VisualStudio::CppUnitTestFramework::Logger::WriteMessage(L"\n");
}
// Equivalent to Log(message); Newline();
void LogLine(std::wstring const& message) {
  WriteLog(message);
  NewLine();
}
// The Microsoft equivalent only takes a wchar_t*.
void AssertTrue(bool const test, std::wstring const& message) {
  Microsoft::VisualStudio::CppUnitTestFramework::Assert::IsTrue(
      test,
      message.c_str());
}

}  // namespace TestUtilities
}  // namespace Principia
