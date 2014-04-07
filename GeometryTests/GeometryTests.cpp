#include "stdafx.hpp"
#include "CppUnitTest.h"
#include "..\Quantities\Quantities.hpp"
#include "..\Quantities\SI.hpp"
#include "..\Geometry\R3Element.hpp"
#include "..\Geometry\Grassmann.hpp"
#include "..\Quantities\ElementaryFunctions.hpp"

using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace Principia {
namespace GeometryTests
{
using namespace Geometry;
using namespace Quantities;
using namespace SI;

TEST_CLASS(GeometryTest)
{
public:
  
  TEST_METHOD(MessyMiscellaneousTests)
  {
    struct World;
    Vector<Length, World> v(R3Element<Length>(0 * Metre, 0 * Metre, 0 * Metre));
    Vector<Length, World> w = -v;
    Length norm = Sqrt(InnerProduct(v, v));
    Bivector<Area, World> s = Wedge(v, w);
  }

};
}

}
