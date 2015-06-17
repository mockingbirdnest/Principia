#pragma once

#include "geometry/serialization.hpp"

#include "base/not_null.hpp"
//#include "geometry/grassmann.hpp"
#include "quantities/quantities.hpp"

namespace principia {

using base::not_null;
using quantities::Quantity;

namespace geometry {

template<typename Message>
class DoubleOrQuantitySerializer<double, Message> {
 public:
  static void WriteToMessage(double const d,
                             not_null<Message*> const message) {
    message->set_double_(d);
  }

  static double ReadFromMessage(Message const& message) {
    CHECK(message.has_double_());
    return message.double_();
  }
};

template<typename Dimensions, typename Message>
class DoubleOrQuantitySerializer<Quantity<Dimensions>, Message> {
 public:
  using T = Quantity<Dimensions>;
  static void WriteToMessage(T const& t, not_null<Message*> const message) {
    t.WriteToMessage(message->mutable_quantity());
  }

  static T ReadFromMessage(Message const& message) {
    CHECK(message.has_quantity());
    return T::ReadFromMessage(message.quantity());
  }
};

template<typename Message>
class DoubleOrQuantityOrMultivectorSerializer<double, Message> {
 public:
  static void WriteToMessage(double const d,
                             not_null<Message*> const message) {
    message->set_double_(d);
  }

  static double ReadFromMessage(Message const& message) {
    CHECK(message.has_double_());
    return message.double_();
  }
};

template<typename Dimensions, typename Message>
class DoubleOrQuantityOrMultivectorSerializer<Quantity<Dimensions>, Message> {
 public:
  using T = Quantity<Dimensions>;
  static void WriteToMessage(T const& t, not_null<Message*> const message) {
    t.WriteToMessage(message->mutable_quantity());
  }

  static T ReadFromMessage(Message const& message) {
    CHECK(message.has_quantity());
    return T::ReadFromMessage(message.quantity());
  }
};

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

}  // namespace geometry
}  // namespace principia
