#include "physics/checkpointer.hpp"

#include "base/status_utilities.hpp"
#include "geometry/named_quantities.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "testing_utilities/matchers.hpp"

namespace principia {
namespace physics {

using ::testing::ElementsAre;
using ::testing::Field;
using ::testing::InSequence;
using ::testing::IsEmpty;
using ::testing::MockFunction;
using ::testing::Ref;
using ::testing::Return;
using ::testing::_;
using namespace principia::base::_not_null;
using namespace principia::geometry::_named_quantities;
using namespace principia::quantities::_si;
using namespace principia::testing_utilities::_matchers;

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

  MockFunction<absl::Status(Message::Checkpoint const&)> reader_;
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
  EXPECT_EQ(t1, checkpointer_.newest_checkpoint());
  EXPECT_THAT(checkpointer_.all_checkpoints(), ElementsAre(t1));

  Instant const t2 = t1 + 8 * Second;
  EXPECT_CALL(writer_, Call(_)).Times(0);
  EXPECT_FALSE(checkpointer_.WriteToCheckpointIfNeeded(
      t2,
      /*max_time_between_checkpoints=*/10 * Second));
  EXPECT_EQ(t1, checkpointer_.oldest_checkpoint());
  EXPECT_EQ(t1, checkpointer_.newest_checkpoint());
  EXPECT_THAT(checkpointer_.all_checkpoints(), ElementsAre(t1));

  EXPECT_CALL(writer_, Call(_));
  Instant const t3 = t2 + 3 * Second;
  EXPECT_TRUE(checkpointer_.WriteToCheckpointIfNeeded(
      t3,
      /*max_time_between_checkpoints=*/10 * Second));
  EXPECT_EQ(t1, checkpointer_.oldest_checkpoint());
  EXPECT_EQ(t3, checkpointer_.newest_checkpoint());
  EXPECT_THAT(checkpointer_.all_checkpoints(), ElementsAre(t1, t3));
}

TEST_F(CheckpointerTest, ReadFromOldestCheckpoint) {
  EXPECT_THAT(checkpointer_.ReadFromOldestCheckpoint(),
              StatusIs(absl::StatusCode::kNotFound));

  Instant const t1 = Instant() + 10 * Second;
  EXPECT_CALL(writer_, Call(_));
  checkpointer_.WriteToCheckpoint(t1);

  EXPECT_CALL(reader_, Call(_))
      .WillOnce(Return(absl::CancelledError()))
      .WillOnce(Return(absl::OkStatus()));
  EXPECT_THAT(checkpointer_.ReadFromOldestCheckpoint(),
              StatusIs(absl::StatusCode::kCancelled));
  EXPECT_OK(checkpointer_.ReadFromOldestCheckpoint());
}

TEST_F(CheckpointerTest, ReadFromNewestCheckpoint) {
  EXPECT_THAT(checkpointer_.ReadFromNewestCheckpoint(),
              StatusIs(absl::StatusCode::kNotFound));

  Instant const t1 = Instant() + 10 * Second;
  EXPECT_CALL(writer_, Call(_));
  checkpointer_.WriteToCheckpoint(t1);

  EXPECT_CALL(reader_, Call(_))
      .WillOnce(Return(absl::CancelledError()))
      .WillOnce(Return(absl::OkStatus()));
  EXPECT_THAT(checkpointer_.ReadFromNewestCheckpoint(),
              StatusIs(absl::StatusCode::kCancelled));
  EXPECT_OK(checkpointer_.ReadFromNewestCheckpoint());
}

TEST_F(CheckpointerTest, ReadFromCheckpointAtOrBefore) {
  Instant const t1 = Instant() + 10 * Second;
  EXPECT_CALL(writer_, Call(_)).WillOnce(SetPayload(1));
  checkpointer_.WriteToCheckpoint(t1);

  Instant const t2 = t1 + 11 * Second;
  EXPECT_CALL(writer_, Call(_)).WillOnce(SetPayload(2));
  checkpointer_.WriteToCheckpoint(t2);

  Instant const t3 = t2 + 11 * Second;
  EXPECT_CALL(writer_, Call(_)).WillOnce(SetPayload(3));
  checkpointer_.WriteToCheckpoint(t3);

  EXPECT_EQ(t1, checkpointer_.checkpoint_at_or_after(t1));
  EXPECT_EQ(t3, checkpointer_.checkpoint_at_or_after(t2 + 1 * Second));
  EXPECT_EQ(InfiniteFuture,
            checkpointer_.checkpoint_at_or_after(t3 + 1 * Second));

  EXPECT_EQ(InfinitePast,
            checkpointer_.checkpoint_at_or_before(Instant() + 1 * Second));
  EXPECT_EQ(t1, checkpointer_.checkpoint_at_or_before(t1));
  EXPECT_EQ(t2, checkpointer_.checkpoint_at_or_before(t2 + 1 * Second));

  EXPECT_THAT(
      checkpointer_.all_checkpoints_at_or_before(Instant() + 1 * Second),
      IsEmpty());
  EXPECT_THAT(checkpointer_.all_checkpoints_at_or_before(t1),
              ElementsAre(t1));
  EXPECT_THAT(checkpointer_.all_checkpoints_at_or_before(t2 + 1 * Second),
              ElementsAre(t1, t2));

  EXPECT_THAT(checkpointer_.all_checkpoints_between(Instant() + 1 * Second,
                                                    Instant() + 3 * Second),
              IsEmpty());
  EXPECT_THAT(checkpointer_.all_checkpoints_between(Instant() + 1 * Second,
                                                    t1),
              ElementsAre(t1));
  EXPECT_THAT(checkpointer_.all_checkpoints_between(Instant() + 1 * Second,
                                                    t2 + 1 * Second),
              ElementsAre(t1, t2));
  EXPECT_THAT(checkpointer_.all_checkpoints_between(t1, t2),
              ElementsAre(t1, t2));
  EXPECT_THAT(checkpointer_.all_checkpoints_between(t1 - 1 * Second,
                                                    t2 + 1 * Second),
              ElementsAre(t1, t2));
  EXPECT_THAT(checkpointer_.all_checkpoints_between(t3, t1),
              IsEmpty());
  EXPECT_THAT(checkpointer_.all_checkpoints_between(t2, t2),
              ElementsAre(t2));
  EXPECT_THAT(checkpointer_.all_checkpoints_between(t1 + 1 * Second,
                                                    t1 + 1 * Second),
              IsEmpty());

  EXPECT_THAT(
      checkpointer_.ReadFromCheckpointAtOrBefore(Instant() + 1 * Second),
      StatusIs(absl::StatusCode::kNotFound));

  EXPECT_CALL(reader_, Call(Field(&Message::Checkpoint::payload, 1)));
  EXPECT_OK(checkpointer_.ReadFromCheckpointAtOrBefore(t1));

  EXPECT_CALL(reader_, Call(Field(&Message::Checkpoint::payload, 2)));
  EXPECT_OK(checkpointer_.ReadFromCheckpointAtOrBefore(t2 + 1 * Second));

  EXPECT_CALL(reader_, Call(Field(&Message::Checkpoint::payload, 3)))
      .WillOnce(Return(absl::CancelledError()));
  EXPECT_THAT(checkpointer_.ReadFromCheckpointAtOrBefore(t3),
              StatusIs(absl::StatusCode::kCancelled));
}

TEST_F(CheckpointerTest, ReadFromCheckpointAt) {
  Instant const t1 = Instant() + 10 * Second;
  EXPECT_CALL(writer_, Call(_)).WillOnce(SetPayload(1));
  checkpointer_.WriteToCheckpoint(t1);

  Instant const t2 = t1 + 11 * Second;
  EXPECT_CALL(writer_, Call(_)).WillOnce(SetPayload(2));
  checkpointer_.WriteToCheckpoint(t2);

  Instant const t3 = t2 + 11 * Second;
  EXPECT_CALL(writer_, Call(_)).WillOnce(SetPayload(3));
  checkpointer_.WriteToCheckpoint(t3);

  EXPECT_THAT(checkpointer_.ReadFromCheckpointAt(t1 + 1 * Second,
                                                 reader_.AsStdFunction()),
              StatusIs(absl::StatusCode::kNotFound));

  EXPECT_CALL(reader_, Call(Field(&Message::Checkpoint::payload, 2)));
  EXPECT_OK(checkpointer_.ReadFromCheckpointAt(t2, reader_.AsStdFunction()));

  EXPECT_CALL(reader_, Call(Field(&Message::Checkpoint::payload, 3)))
      .WillOnce(Return(absl::CancelledError()));
  EXPECT_THAT(checkpointer_.ReadFromCheckpointAt(t3, reader_.AsStdFunction()),
              StatusIs(absl::StatusCode::kCancelled));
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
