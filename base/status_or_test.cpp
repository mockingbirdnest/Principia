
// This code comes from:
// https://github.com/google/protobuf/tree/master/src/google/protobuf/stubs
// and was adapted to Visual Studio and to the needs of this project.

// Protocol Buffers - Google's data interchange format
// Copyright 2008 Google Inc.  All rights reserved.
// https://developers.google.com/protocol-buffers/
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above
// copyright notice, this list of conditions and the following disclaimer
// in the documentation and/or other materials provided with the
// distribution.
//     * Neither the name of Google Inc. nor the names of its
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include "base/status_or.hpp"

#include <cerrno>
#include <memory>

#include "gtest/gtest.h"

namespace principia {
namespace base {

class Base1 {
 public:
  virtual ~Base1() = default;
 private:
  [[maybe_unused]] int pad;
};

class Base2 {
 public:
  virtual ~Base2() = default;
 private:
  [[maybe_unused]] int yetotherpad;
};

class Derived : public Base1, public Base2 {
 public:
  ~Derived() override = default;
 private:
  [[maybe_unused]] int evenmorepad;
};

TEST(StatusOr, TestDefaultCtor) {
  StatusOr<int> thing;
  EXPECT_FALSE(thing.ok());
  EXPECT_EQ(Status::UNKNOWN, thing.status());
}

TEST(StatusOr, TestStatusCtor) {
  StatusOr<int> thing(Status::CANCELLED);
  EXPECT_FALSE(thing.ok());
  EXPECT_EQ(Status::CANCELLED, thing.status());
}

TEST(StatusOr, TestValueCtor) {
  const int i = 4;
  StatusOr<int> thing(i);
  EXPECT_TRUE(thing.ok());
  EXPECT_EQ(i, thing.ValueOrDie());
}

TEST(StatusOr, TestCopyCtorStatusOk) {
  const int i = 4;
  StatusOr<int> original(i);
  StatusOr<int> copy(original);
  EXPECT_EQ(original.status(), copy.status());
  EXPECT_EQ(original.ValueOrDie(), copy.ValueOrDie());
}

TEST(StatusOr, TestCopyCtorStatusNotOk) {
  StatusOr<int> original(Status::CANCELLED);
  StatusOr<int> copy(original);
  EXPECT_EQ(original.status(), copy.status());
}

TEST(StatusOr, TestCopyCtorStatusOKConverting) {
  const int i = 4;
  StatusOr<int> original(i);
  StatusOr<double> copy(original);
  EXPECT_EQ(original.status(), copy.status());
  EXPECT_EQ(original.ValueOrDie(), copy.ValueOrDie());
}

TEST(StatusOr, TestCopyCtorStatusNotOkConverting) {
  StatusOr<int> original(Status::CANCELLED);
  StatusOr<double> copy(original);
  EXPECT_EQ(original.status(), copy.status());
}

TEST(StatusOr, TestAssignmentStatusOk) {
  const int i = 4;
  StatusOr<int> source(i);
  StatusOr<int> target;
  target = source;
  EXPECT_EQ(source.status(), target.status());
  EXPECT_EQ(source.ValueOrDie(), target.ValueOrDie());
}

TEST(StatusOr, TestAssignmentStatusNotOk) {
  StatusOr<int> source(Status::CANCELLED);
  StatusOr<int> target;
  target = source;
  EXPECT_EQ(source.status(), target.status());
}

TEST(StatusOr, TestAssignmentStatusOKConverting) {
  const int i = 4;
  StatusOr<int> source(i);
  StatusOr<double> target;
  target = source;
  EXPECT_EQ(source.status(), target.status());
  EXPECT_DOUBLE_EQ(source.ValueOrDie(), target.ValueOrDie());
}

TEST(StatusOr, TestAssignmentStatusNotOkConverting) {
  StatusOr<int> source(Status::CANCELLED);
  StatusOr<double> target;
  target = source;
  EXPECT_EQ(source.status(), target.status());
}

TEST(StatusOr, TestStatus) {
  StatusOr<int> good(4);
  EXPECT_TRUE(good.ok());
  StatusOr<int> bad(Status::CANCELLED);
  EXPECT_FALSE(bad.ok());
  EXPECT_EQ(Status::CANCELLED, bad.status());
}

TEST(StatusOr, TestValue) {
  const int i = 4;
  StatusOr<int> thing(i);
  EXPECT_EQ(i, thing.ValueOrDie());
}

TEST(StatusOr, TestValueConst) {
  const int i = 4;
  const StatusOr<int> thing(i);
  EXPECT_EQ(i, thing.ValueOrDie());
}

TEST(StatusOr, TestPointerDefaultCtor) {
  StatusOr<int*> thing;
  EXPECT_FALSE(thing.ok());
  EXPECT_EQ(Status::UNKNOWN, thing.status());
}

TEST(StatusOr, TestPointerStatusCtor) {
  StatusOr<int*> thing(Status::CANCELLED);
  EXPECT_FALSE(thing.ok());
  EXPECT_EQ(Status::CANCELLED, thing.status());
}

TEST(StatusOr, TestPointerValueCtor) {
  const int i = 4;
  StatusOr<const int*> thing(&i);
  EXPECT_TRUE(thing.ok());
  EXPECT_EQ(&i, thing.ValueOrDie());
}

TEST(StatusOr, TestPointerCopyCtorStatusOk) {
  const int i = 0;
  StatusOr<const int*> original(&i);
  StatusOr<const int*> copy(original);
  EXPECT_EQ(original.status(), copy.status());
  EXPECT_EQ(original.ValueOrDie(), copy.ValueOrDie());
}

TEST(StatusOr, TestPointerCopyCtorStatusNotOk) {
  StatusOr<int*> original(Status::CANCELLED);
  StatusOr<int*> copy(original);
  EXPECT_EQ(original.status(), copy.status());
}

TEST(StatusOr, TestPointerCopyCtorStatusOKConverting) {
  Derived derived;
  StatusOr<Derived*> original(&derived);
  StatusOr<Base2*> copy(original);
  EXPECT_EQ(original.status(), copy.status());
  EXPECT_EQ(static_cast<const Base2*>(original.ValueOrDie()),
            copy.ValueOrDie());
}

TEST(StatusOr, TestPointerCopyCtorStatusNotOkConverting) {
  StatusOr<Derived*> original(Status::CANCELLED);
  StatusOr<Base2*> copy(original);  // NOLINT
  EXPECT_EQ(original.status(), copy.status());
}

TEST(StatusOr, TestPointerAssignmentStatusOk) {
  const int i = 0;
  StatusOr<const int*> source(&i);
  StatusOr<const int*> target;
  target = source;
  EXPECT_EQ(source.status(), target.status());
  EXPECT_EQ(source.ValueOrDie(), target.ValueOrDie());
}

TEST(StatusOr, TestPointerAssignmentStatusNotOk) {
  StatusOr<int*> source(Status::CANCELLED);
  StatusOr<int*> target;
  target = source;
  EXPECT_EQ(source.status(), target.status());
}

TEST(StatusOr, TestPointerAssignmentStatusOKConverting) {
  Derived derived;
  StatusOr<Derived*> source(&derived);
  StatusOr<Base2*> target;
  target = source;
  EXPECT_EQ(source.status(), target.status());
  EXPECT_EQ(static_cast<const Base2*>(source.ValueOrDie()),
            target.ValueOrDie());
}

TEST(StatusOr, TestPointerAssignmentStatusNotOkConverting) {
  StatusOr<Derived*> source(Status::CANCELLED);
  StatusOr<Base2*> target;
  target = source;
  EXPECT_EQ(source.status(), target.status());
}

TEST(StatusOr, TestPointerStatus) {
  const int i = 0;
  StatusOr<const int*> good(&i);
  EXPECT_TRUE(good.ok());
  StatusOr<const int*> bad(Status::CANCELLED);
  EXPECT_EQ(Status::CANCELLED, bad.status());
}

TEST(StatusOr, TestPointerValue) {
  const int i = 0;
  StatusOr<const int*> thing(&i);
  EXPECT_EQ(&i, thing.ValueOrDie());
}

TEST(StatusOr, TestPointerValueConst) {
  const int i = 0;
  const StatusOr<const int*> thing(&i);
  EXPECT_EQ(&i, thing.ValueOrDie());
}

}  // namespace base
}  // namespace principia
