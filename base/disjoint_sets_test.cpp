#include "base/disjoint_sets.hpp"

#include <vector>

#include "gtest/gtest.h"
#include "gmock/gmock.h"

using ::testing::Eq;

namespace principia {
namespace base {

using namespace principia::base::_disjoint_sets;

namespace {

struct CountableInteger;

}  // namespace

template<>
class Subset<CountableInteger>::Properties final {
 public:
  explicit Properties(int const cardinality) : cardinality(cardinality) {}

  void MergeWith(Properties& other) {
    cardinality += other.cardinality;
    other.cardinality = -1;
  }

  int cardinality = -1;
};

namespace {

struct Integer final {
  int value;
  Subset<Integer>::Node subset_node;
};

struct CountableInteger final {
  int value;
  Subset<CountableInteger>::Node subset_node;
};

}  // namespace

class DisjointSetsTest : public testing::Test {
 protected:
  DisjointSetsTest() {
    integers_.resize(50);
    countable_integers_.resize(50);
    for (int i = 0; i < 50; ++i) {
      integers_[i].value = i;
      countable_integers_[i].value = i;
      Subset<Integer>::MakeSingleton(integers_[i]);
      Subset<CountableInteger>::MakeSingleton(countable_integers_[i], 1);
    }
  }

  std::vector<Integer> integers_;
  std::vector<CountableInteger> countable_integers_;
};


template<>
not_null<Subset<Integer>::Node*> Subset<Integer>::Node::Get(Integer& element) {
  return &element.subset_node;
}

template<>
not_null<Subset<CountableInteger>::Node*> Subset<CountableInteger>::Node::Get(
    CountableInteger& element) {
  return &element.subset_node;
}

TEST_F(DisjointSetsTest, Congruence) {
  for (auto& left : integers_) {
    for (auto& right : integers_) {
      if (left.value % 5 == right.value % 5) {
        auto const unified =
            Subset<Integer>::Unite(Subset<Integer>::Find(left),
                                   Subset<Integer>::Find(right));
        EXPECT_EQ(unified, Subset<Integer>::Find(left));
        EXPECT_EQ(unified, Subset<Integer>::Find(right));
      }
    }
  }
  for (auto& i : integers_) {
    EXPECT_EQ(Subset<Integer>::Find(integers_[i.value % 5]),
              Subset<Integer>::Find(i));
  }

  for (auto& left : countable_integers_) {
    for (auto& right : countable_integers_) {
      if (left.value % 5 == right.value % 5) {
        Subset<CountableInteger>::Unite(Subset<CountableInteger>::Find(left),
                                        Subset<CountableInteger>::Find(right));
      }
    }
  }
  for (auto& i : countable_integers_) {
    EXPECT_EQ(10, Subset<CountableInteger>::Find(i).properties().cardinality);
  }
}

}  // namespace base
}  // namespace principia
