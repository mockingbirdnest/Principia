#include "base/disjoint_sets.hpp"

#include <vector>

#include "gtest/gtest.h"
#include "gmock/gmock.h"

using testing::Eq;

namespace principia {
namespace base {

namespace {

struct CountableInteger;

}  // namespace

template<>
class SubsetProperties<CountableInteger> {
 public:
  explicit SubsetProperties<CountableInteger>(int const cardinality)
      : cardinality(cardinality) {}

  void MergeWith(SubsetProperties& other) {
    cardinality += other.cardinality;
    other.cardinality = -1;
  }

  int cardinality = -1;
};

namespace {

struct Integer {
  int value;
  SubsetNode<Integer> subset_node;
};

struct CountableInteger {
  int value;
  SubsetNode<CountableInteger> subset_node;
};

}  // namespace

class DisjointSetsTest : public testing::Test {
 protected:
  DisjointSetsTest() {
    integers_.resize(50);
    countable_integers_.resize(50);
    for(int i = 0; i < 50; ++i) {
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
not_null<SubsetNode<Integer>*> GetSubsetNode<Integer>(Integer& element) {
  return &element.subset_node;
}

template<>
not_null<SubsetNode<CountableInteger>*> GetSubsetNode<CountableInteger>(
    CountableInteger& element) {
  return &element.subset_node;
}

TEST_F(DisjointSetsTest, Congruence) {
  for (auto& left : integers_) {
    for (auto& right : integers_) {
      if (left.value % 5 == right.value % 5) {
        auto const unified =
            Unite(Subset<Integer>::Find(left), Subset<Integer>::Find(right));
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
        Unite(Subset<CountableInteger>::Find(left),
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
