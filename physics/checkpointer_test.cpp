#include "physics/checkpointer.hpp"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

namespace principia {
namespace physics {

using base::not_null;
using geometry::Instant;
using quantities::si::Second;
using ::testing::MockFunction;
using ::testing::Ref;
using ::testing::Return;
using ::testing::_;

struct Message {
  Message() {}
  Message(const Message&) {}
  MOCK_METHOD1(MergeFrom, void(Message const&));
};

class CheckpointerTest : public ::testing::Test {
 protected:
  CheckpointerTest()
      : checkpointer_(reader_.AsStdFunction(),
                      writer_.AsStdFunction()) {}

  MockFunction<bool(Message const&)> reader_;
  MockFunction<void(not_null<Message*>)> writer_;
  Checkpointer<Message> checkpointer_;
};

TEST_F(CheckpointerTest, CreateUnconditionally) {
  Instant const t = Instant() + 10 * Second;
  EXPECT_CALL(writer_, Call(_));
  checkpointer_.CreateUnconditionally(t);
  Message m;
  EXPECT_CALL(m, MergeFrom(_));
  EXPECT_EQ(t, checkpointer_.WriteToMessage(&m));
}

TEST_F(CheckpointerTest, CreateIfNeeded) {
  Instant const t1 = Instant() + 10 * Second;
  EXPECT_CALL(writer_, Call(_));
  checkpointer_.CreateUnconditionally(t1);

  Instant const t2 = t1 + 8 * Second;
  EXPECT_CALL(writer_, Call(_)).Times(0);
  checkpointer_.CreateIfNeeded(t2,
                               /*max_time_between_checkpoints=*/10 * Second);

  EXPECT_CALL(writer_, Call(_));
  Instant const t3 = t2 + 3 * Second;
  checkpointer_.CreateIfNeeded(t3,
                               /*max_time_between_checkpoints=*/10 * Second);

  Message m;
  EXPECT_CALL(m, MergeFrom(_));
  EXPECT_EQ(t1, checkpointer_.WriteToMessage(&m));
}

TEST_F(CheckpointerTest, ForgetBefore) {
  Instant const t1 = Instant() + 10 * Second;
  EXPECT_CALL(writer_, Call(_));
  checkpointer_.CreateUnconditionally(t1);

  Instant const t2 = t1 + 8 * Second;
  EXPECT_CALL(writer_, Call(_));
  checkpointer_.CreateUnconditionally(t2);

  checkpointer_.ForgetBefore(t1 + 4 * Second);

  Message m;
  EXPECT_CALL(m, MergeFrom(_));
  EXPECT_EQ(t2, checkpointer_.WriteToMessage(&m));
}

TEST_F(CheckpointerTest, ReadFromMessage) {
  Instant const t = Instant() + 10 * Second;
  Message m;

  EXPECT_CALL(reader_, Call(Ref(m))).WillOnce(Return(false));
  checkpointer_.ReadFromMessage(t, m);

  EXPECT_CALL(reader_, Call(Ref(m))).WillOnce(Return(true));
  EXPECT_CALL(writer_, Call(_));
  checkpointer_.ReadFromMessage(t, m);
}

TEST_F(CheckpointerTest, WriteToMessage) {
  Instant const t = Instant() + 10 * Second;
  Message m;

  EXPECT_CALL(m, MergeFrom(_)).Times(0);
  EXPECT_LT(t + 1000 * Second, checkpointer_.WriteToMessage(&m));

  EXPECT_CALL(writer_, Call(_));
  checkpointer_.CreateUnconditionally(t);

  EXPECT_CALL(m, MergeFrom(_));
  EXPECT_EQ(t, checkpointer_.WriteToMessage(&m));
}

}  // namespace physics
}  // namespace principia
