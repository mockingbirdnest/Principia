
#include "geometry/symmetric_bilinear_form.hpp"

#include "geometry/frame.hpp"
#include "gtest/gtest.h"
#include "quantities/quantities.hpp"
#include "serialization/geometry.pb.h"

namespace principia {
namespace geometry {
namespace internal_symmetric_bilinear_form {

using quantities::Length;

class SymmetricBilinearFormTest : public ::testing::Test {
 protected:
  using World =
      Frame<serialization::Frame::TestTag, serialization::Frame::TEST, true>;

  SymmetricBilinearFormTest()
      : form_(SymmetricBilinearForm<Length, World>::InnerProductForm()) {}

  SymmetricBilinearForm<Length, World> form_;
};

}  // namespace internal_symmetric_bilinear_form
}  // namespace geometry
}  // namespace principia
