#pragma once

#include "geometry/serialization.hpp"

#include "base/not_null.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/point.hpp"
#include "quantities/quantities.hpp"
#include "quantities/serialization.hpp"

namespace principia {

using base::not_null;
using quantities::DoubleOrQuantitySerializer;
using quantities::Quantity;

namespace geometry {

template<typename Message>
class DoubleOrQuantityOrMultivectorSerializer<double, Message>
    : public DoubleOrQuantitySerializer<double, Message> {};

template<typename Dimensions, typename Message>
class DoubleOrQuantityOrMultivectorSerializer<Quantity<Dimensions>, Message>
    : public DoubleOrQuantitySerializer<Quantity<Dimensions>, Message> {};

template<typename Scalar, typename Frame, int rank, typename Message>
class DoubleOrQuantityOrMultivectorSerializer<
    Multivector<Scalar, Frame, rank>, Message> {
 public:
  using T = Multivector<Scalar, Frame, rank>;
  static void WriteToMessage(T const& t, not_null<Message*> const message) {
    t.WriteToMessage(message->mutable_multivector());
  }

  static T ReadFromMessage(Message const& message) {
    CHECK(message.has_multivector());
    return T::ReadFromMessage(message.multivector());
  }
};

template<typename Vector, typename Message>
class PointOrMultivectorSerializer<Point<Vector>, Message> {
 public:
  using T = Point<Vector>;
  static void WriteToMessage(T const& t, not_null<Message*> const message) {
    t.WriteToMessage(message->mutable_point());
  }

  static T ReadFromMessage(Message const& message) {
    CHECK(message.has_point());
    return T::ReadFromMessage(message.point());
  }
};

template<typename Scalar, typename Frame, int rank, typename Message>
class PointOrMultivectorSerializer<Multivector<Scalar, Frame, rank>, Message>
    : public DoubleOrQuantityOrMultivectorSerializer<
                 Multivector<Scalar, Frame, rank>, Message> {};

template<typename Scalar, typename Frame, int rank, typename Message>
class QuantityOrMultivectorSerializer<Multivector<Scalar, Frame, rank>, Message>
    : public DoubleOrQuantityOrMultivectorSerializer<
                 Multivector<Scalar, Frame, rank>, Message> {};

template<typename Dimensions, typename Message>
class QuantityOrMultivectorSerializer<Quantity<Dimensions>, Message>
    : public DoubleOrQuantitySerializer<Quantity<Dimensions>, Message> {};

}  // namespace geometry
}  // namespace principia
