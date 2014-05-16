#pragma once

#include <CppUnitTest.h>

#include "quantities/dimensionless.hpp"
#include "quantities/quantities.hpp"

namespace principia {
namespace test_utilities {

template<typename ValueType>
void AssertEqualWithin(ValueType const& left,
                       ValueType const& right,
                       quantities::Dimensionless const& ε) {
  std::string message = "Should be equal within " + ToString(ε, 3) +
                         ": " + ToString(left) + " and " + ToString(right) +
                         ".";
  LogLine(message);
  AssertTrue(left == right || Abs(left / right - 1) < ε, message);
}

template<typename ValueType>
void AssertEqual(ValueType const& left,
                 ValueType const& right) {
  std::string message = "Should be equal: " + ToString(left) + " and " +
                        ToString(right) + ".";
  LogLine(message);
  AssertTrue(left == right, message);
}

template<typename ValueType>
void AssertNotEqualWithin(ValueType const& left,
                          ValueType const& right,
                          quantities::Dimensionless const& ε) {
  std::string message = "Should differ by more than " + ToString(ε, 3) +
                         ": " + ToString(left) + " and " + ToString(right) +
                         ".";
  LogLine(message);
  AssertTrue(Abs(left / right - 1) > ε, message);
}

template<typename ValueType>
void AssertNotEqual(ValueType const& left, ValueType const& right) {
  std::string message = "Should be different: " + ToString(left) + " and " +
                        ToString(right) + ".";
  LogLine(message);
  AssertTrue(left != right, message);
}

inline void WriteLog(std::string const& message) {
  Microsoft::VisualStudio::CppUnitTestFramework::Logger::WriteMessage(
      message.c_str());
}

inline void NewLine() {
  Microsoft::VisualStudio::CppUnitTestFramework::Logger::WriteMessage("\n");
}

inline void LogLine(std::string const& message) {
  WriteLog(message);
  NewLine();
}

inline void AssertFalse(bool const test, std::string const& message) {
  Microsoft::VisualStudio::CppUnitTestFramework::Assert::IsFalse(
      test,
      std::wstring(message.begin(), message.end()).c_str());
}

inline void AssertTrue(bool const test, std::string const& message) {
  Microsoft::VisualStudio::CppUnitTestFramework::Assert::IsTrue(
      test,
      std::wstring(message.begin(), message.end()).c_str());
}

}  // namespace test_utilities
}  // namespace principia
