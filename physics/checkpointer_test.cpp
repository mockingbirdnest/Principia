#include "physics/checkpointer.hpp"

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "testing_utilities/matchers.hpp"

namespace principia {
namespace physics {

using base::Error;
using base::not_null;
using base::Status;
using geometry::Instant;
using quantities::si::Second;
using testing_utilities::StatusIs;
using ::testing::Field;
using ::testing::InSequence;
using ::testing::MockFunction;
using ::testing::Ref;
using ::testing::Return;
using ::testing::_;

ACTION_P(SetPayload, payload) {
  arg0->payload = payload;
}

struct Message {
  class Checkpoint {
   public:
    serialization::Point* mutable_time() {
      return &time_;
    }
    const serialization::Point& time() const {
      return time_;
    }

    int payload = 0;

   private:
    serialization::Point time_;
  };
  google::protobuf::RepeatedPtrField<Checkpoint> checkpoint;
};

class CheckpointerTest : public ::testing::Test {
 protected:
  CheckpointerTest()
      : checkpointer_(writer_.AsStdFunction(),
                      reader_.AsStdFunction()) {}

  MockFunction<Status(Message::Checkpoint const&)> reader_;
  MockFunction<void(not_null<Message::Checkpoint*>)> writer_;
  Checkpointer<Message> checkpointer_;
};

TEST_F(CheckpointerTest, WriteToCheckpoint) {
  Instant const t = Instant() + 10 * Second;
  EXPECT_CALL(writer_, Call(_));
  checkpointer_.WriteToCheckpoint(t);
}

TEST_F(CheckpointerTest, WriteToCheckpointIfNeeded) {
  Instant const t1 = Instant() + 10 * Second;
  EXPECT_CALL(writer_, Call(_));
  checkpointer_.WriteToCheckpoint(t1);
  EXPECT_EQ(t1, checkpointer_.oldest_checkpoint());

  Instant const t2 = t1 + 8 * Second;
  EXPECT_CALL(writer_, Call(_)).Times(0);
  EXPECT_FALSE(checkpointer_.WriteToCheckpointIfNeeded(
      t2,
      /*max_time_between_checkpoints=*/10 * Second));
  EXPECT_EQ(t1, checkpointer_.oldest_checkpoint());

  EXPECT_CALL(writer_, Call(_));
  Instant const t3 = t2 + 3 * Second;
  EXPECT_TRUE(checkpointer_.WriteToCheckpointIfNeeded(
      t3,
      /*max_time_between_checkpoints=*/10 * Second));
  EXPECT_EQ(t1, checkpointer_.oldest_checkpoint());
}

TEST_F(CheckpointerTest, ReadFromOldestCheckpoint) {
  EXPECT_THAT(checkpointer_.ReadFromOldestCheckpoint(),
              StatusIs(Error::NOT_FOUND));

  Instant const t1 = Instant() + 10 * Second;
  EXPECT_CALL(writer_, Call(_));
  checkpointer_.WriteToCheckpoint(t1);

  EXPECT_CALL(reader_, Call(_))
      .WillOnce(Return(Status::CANCELLED))
      .WillOnce(Return(Status::OK));
  EXPECT_THAT(checkpointer_.ReadFromOldestCheckpoint(),
              StatusIs(Error::CANCELLED));
  EXPECT_OK(checkpointer_.ReadFromOldestCheckpoint());
}

TEST_F(CheckpointerTest, ReadFromAllCheckpointsBackwards) {
  Instant const t1 = Instant() + 10 * Second;
  EXPECT_CALL(writer_, Call(_)).WillOnce(SetPayload(1));
  checkpointer_.WriteToCheckpoint(t1);

  Instant const t2 = t1 + 11 * Second;
  EXPECT_CALL(writer_, Call(_)).WillOnce(SetPayload(2));
  checkpointer_.WriteToCheckpoint(t2);

  Instant const t3 = t2 + 11 * Second;
  EXPECT_CALL(writer_, Call(_)).WillOnce(SetPayload(3));
  checkpointer_.WriteToCheckpoint(t3);

  {
    InSequence s;
    EXPECT_CALL(reader_, Call(Field(&Message::Checkpoint::payload, 3)));
    EXPECT_CALL(reader_, Call(Field(&Message::Checkpoint::payload, 2)));
    EXPECT_CALL(reader_, Call(Field(&Message::Checkpoint::payload, 1)));
  }
  EXPECT_OK(
      checkpointer_.ReadFromAllCheckpointsBackwards(reader_.AsStdFunction()));

  {
    InSequence s;
    EXPECT_CALL(reader_, Call(Field(&Message::Checkpoint::payload, 3)));
    EXPECT_CALL(reader_, Call(Field(&Message::Checkpoint::payload, 2)))
        .WillOnce(Return(Status::CANCELLED));
  }
  EXPECT_THAT(
      checkpointer_.ReadFromAllCheckpointsBackwards(reader_.AsStdFunction()),
      StatusIs(Error::CANCELLED));
}

TEST_F(CheckpointerTest, Serialization) {
  Instant t = Instant() + 10 * Second;
  EXPECT_CALL(writer_, Call(_)).Times(2);
  checkpointer_.WriteToCheckpoint(t);
  t += 13 * Second;
  checkpointer_.WriteToCheckpoint(t);

  Message m;
  checkpointer_.WriteToMessage(&m.checkpoint);
  EXPECT_EQ(2, m.checkpoint.size());
  EXPECT_EQ(10, m.checkpoint[0].time().scalar().magnitude());
  EXPECT_EQ(23, m.checkpoint[1].time().scalar().magnitude());

  auto const checkpointer =
      Checkpointer<Message>::ReadFromMessage(writer_.AsStdFunction(),
                                             reader_.AsStdFunction(),
                                             m.checkpoint);
  EXPECT_EQ(Instant() + 10 * Second, checkpointer->oldest_checkpoint());
}

}  // namespace physics
}  // namespace principia
