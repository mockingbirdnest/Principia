#include "stdafx.h"
#include "CppUnitTest.h"

#include "..\PhysicalQuantities\PhysicalQuantities.h"
#include "..\PhysicalQuantities\NamedQuantities.h"
#include "..\PhysicalQuantities\SIUnits.h"
#include "..\PhysicalQuantities\Constants.h"
#include <stdio.h>

using namespace Microsoft::VisualStudio::CppUnitTestFramework;

using namespace PhysicalQuantities;

namespace QuantitiesTests {		
TEST_CLASS(UnitTest1) {
 public:
  TEST_METHOD(PhysicalConstants) {
    Dimensionless ε = 1e-15;
    Assert::IsTrue(Abs((1 / SpeedOfLight.Pow<2>()) / 
                       (ElectricConstant * MagneticConstant) - 1) < ε);
    Logger::WriteMessage(ToString(1 / SpeedOfLight.Pow<2>()).c_str());
    Logger::WriteMessage(L"\n");
    Logger::WriteMessage(ToString(MagneticConstant * ElectricConstant).c_str());
  }
};
}