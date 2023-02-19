#include "physics/clientele.hpp"

#include <memory>

#include "gtest/gtest.h"

namespace principia {
namespace physics {

class ClienteleTest : public ::testing::Test {
 protected:
  ClienteleTest() : clientele_(/*default_key=*/666) {}

  Clientele<int> clientele_;
};

TEST_F(ClienteleTest, Clientele) {
  EXPECT_EQ(666, clientele_.first());
  clientele_.Join(1);
  clientele_.Join(2);
  clientele_.Join(1);
  EXPECT_EQ(1, clientele_.first());
  clientele_.Leave(1);
  EXPECT_EQ(1, clientele_.first());
  clientele_.Leave(1);
  EXPECT_EQ(2, clientele_.first());
}

TEST_F(ClienteleTest, Client) {
  std::unique_ptr<Client<int>> client0;
  {
    Client const client2(2, clientele_);
    EXPECT_EQ(2, clientele_.first());
    {
      Client const client1(1, clientele_);
      EXPECT_EQ(1, clientele_.first());
    }
    EXPECT_EQ(2, clientele_.first());
    client0 = std::make_unique<Client<int>>(0, clientele_);
    EXPECT_EQ(0, clientele_.first());
  }
  EXPECT_EQ(0, clientele_.first());
}

}  // namespace physics
}  // namespace principia
