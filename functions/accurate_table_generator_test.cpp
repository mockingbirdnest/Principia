#include "functions/accurate_table_generator.hpp"

#include "functions/multiprecision.hpp"
#include "gtest/gtest.h"

namespace principia {
namespace functions {
namespace _multiprecision {

using namespace principia::functions::_accurate_table_generator;
using namespace principia::functions::_multiprecision;

class AccurateTableGeneratorTest : public ::testing::Test {};

TEST_F(AccurateTableGeneratorTest, Smoke) {
  ExhaustiveSearch<5>({Sin}, 1.0 / 256.0);
}

}  // namespace _multiprecision
}  // namespace functions
}  // namespace principia
