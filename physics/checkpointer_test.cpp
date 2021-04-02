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
  class Checkpoint {
   public:
    serialization::Point* mutable_time() {
      return &time_;
    }
    const serialization::Point& time() const {
      return time_;
    }

   private:
    serialization::Point time_;
  };
  google::protobuf::RepeatedPtrField<Checkpoint> checkpoint;
};

class CheckpointerTest : public ::testing::Test {
 protected:
  CheckpointerTest()
      : checkpointer_(reader_.AsStdFunction(),
                      writer_.AsStdFunction()) {}

  MockFunction<bool(Message::Checkpoint const&)> reader_;
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

  Instant const t2 = t1 + 8 * Second;
  EXPECT_CALL(writer_, Call(_)).Times(0);
  checkpointer_.WriteToCheckpointIfNeeded(t2,
                               /*max_time_between_checkpoints=*/10 * Second);

  EXPECT_CALL(writer_, Call(_));
  Instant const t3 = t2 + 3 * Second;
  checkpointer_.WriteToCheckpointIfNeeded(t3,
                               /*max_time_between_checkpoints=*/10 * Second);
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
      Checkpointer<Message>::ReadFromMessage(reader_.AsStdFunction(),
                                             writer_.AsStdFunction(),
                                             m.checkpoint);
  EXPECT_EQ(Instant() + 10 * Second, checkpointer->oldest_checkpoint());
}

}  // namespace physics
}  // namespace principia
