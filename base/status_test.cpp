
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
#include "base/status.hpp"

#include <stdio.h>
#include <string>

#include "gtest/gtest.h"

namespace principia {
namespace base {

TEST(Status, Empty) {
  Status status;
  EXPECT_EQ(Error::OK, Status::OK.error());
  EXPECT_EQ("OK", Status::OK.ToString());
}

TEST(Status, GenericCodes) {
  EXPECT_EQ(Error::OK, Status::OK.error());
  EXPECT_EQ(Error::CANCELLED, Status::CANCELLED.error());
  EXPECT_EQ(Error::UNKNOWN, Status::UNKNOWN.error());
}

TEST(Status, ConstructorZero) {
  Status status(Error::OK, "msg");
  EXPECT_TRUE(status.ok());
  EXPECT_EQ("OK", status.ToString());
}

TEST(Status, CheckOKOK) {
  CHECK_OK(Status::OK);
}

TEST(StatusDeathTest, CheckOKCANCELLED) {
  EXPECT_DEATH({
    CHECK_OK(Status::CANCELLED);
  }, "CANCELLED");
}

TEST(Status, Message) {
  Status status(Error::INVALID_ARGUMENT, "");
  EXPECT_FALSE(status.ok());
  EXPECT_EQ("", status.message());
  EXPECT_EQ("INVALID_ARGUMENT", status.ToString());
  status = Status(Error::INVALID_ARGUMENT, "msg");
  EXPECT_FALSE(status.ok());
  EXPECT_EQ("msg", status.message());
  EXPECT_EQ("INVALID_ARGUMENT: msg", status.ToString());
  status = Status(Error::OK, "msg");
  EXPECT_TRUE(status.ok());
  EXPECT_EQ("", status.message());
  EXPECT_EQ("OK", status.ToString());
}

TEST(Status, Copy) {
  Status a(Error::UNKNOWN, "message");
  Status b(a);
  ASSERT_EQ(a.ToString(), b.ToString());
}

TEST(Status, Assign) {
  Status a(Error::UNKNOWN, "message");
  Status b;
  b = a;
  ASSERT_EQ(a.ToString(), b.ToString());
}

TEST(Status, AssignEmpty) {
  Status a(Error::UNKNOWN, "message");
  Status b;
  a = b;
  ASSERT_EQ(std::string("OK"), a.ToString());
  ASSERT_TRUE(b.ok());
  ASSERT_TRUE(a.ok());
}

TEST(Status, EqualsOK) {
  ASSERT_EQ(Status::OK, Status());
}

TEST(Status, EqualsSame) {
  const Status a = Status(Error::CANCELLED, "message");
  const Status b = Status(Error::CANCELLED, "message");
  ASSERT_EQ(a, b);
}

TEST(Status, EqualsCopy) {
  const Status a = Status(Error::CANCELLED, "message");
  const Status b = a;
  ASSERT_EQ(a, b);
}

TEST(Status, EqualsDifferentCode) {
  const Status a = Status(Error::CANCELLED, "message");
  const Status b = Status(Error::UNKNOWN, "message");
  ASSERT_NE(a, b);
}

TEST(Status, EqualsDifferentMessage) {
  const Status a = Status(Error::CANCELLED, "message");
  const Status b = Status(Error::CANCELLED, "another");
  ASSERT_NE(a, b);
}

}  // namespace base
}  // namespace principia
