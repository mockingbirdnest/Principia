#pragma once

namespace Principia {
namespace Quantities {

namespace TypeGenerators {
  template<typename Q>
struct SquareRootGenerator<
    Q, Condition<! (Q::Dimensions::Length & 1 ||
                    Q::Dimensions::Mass & 1 ||
                    Q::Dimensions::Time & 1 ||
                    Q::Dimensions::Current & 1 ||
                    Q::Dimensions::Temperature & 1 ||
                    Q::Dimensions::Amount & 1 ||
                    Q::Dimensions::LuminousIntensity & 1 ||
                    Q::Dimensions::Winding & 1 ||
                    Q::Dimensions::Angle & 1 ||
                    Q::Dimensions::SolidAngle & 1)>> {
  enum {
    Length            = Q::Dimensions::Length / 2,
    Mass              = Q::Dimensions::Mass / 2,
    Time              = Q::Dimensions::Time / 2,
    Current           = Q::Dimensions::Current / 2,
    Temperature       = Q::Dimensions::Temperature / 2,
    Amount            = Q::Dimensions::Amount / 2,
    LuminousIntensity = Q::Dimensions::LuminousIntensity / 2,
    Winding           = Q::Dimensions::Winding / 2,
    Angle             = Q::Dimensions::Angle / 2,
    SolidAngle        = Q::Dimensions::SolidAngle / 2
  };
  typedef Quantity<
      Dimensions<Length, Mass, Time, Current, Temperature, Amount,
                 LuminousIntensity, Winding, Angle, SolidAngle>> ResultType;
};
}

}
}