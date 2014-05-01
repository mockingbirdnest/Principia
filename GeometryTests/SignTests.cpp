#include "stdafx.hpp"
#include "CppUnitTest.h"
#include "..\Geometry\Sign.hpp"
#include "..\TestUtilities\TestUtilities.hpp"

using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace Principia {
namespace Geometry {
namespace {
using namespace TestUtilities;

TEST_CLASS(SignTests) {
 public:
	TEST_METHOD(Basic) {
    const Sign positive(true);
    const Sign negative(false);
    AssertTrue(positive.Positive());
    AssertFalse(positive.Negative());
    AssertFalse(negative.Positive());
    AssertTrue(negative.Negative());
	}

	TEST_METHOD(Multiplication) {
    const Sign positive(true);
    const Sign negative(false);
    AssertTrue((positive * positive).Positive());
    AssertTrue((positive * negative).Negative());
    AssertTrue((negative * positive).Negative());
    AssertTrue((negative * negative).Positive());
	}
};

}  // namespace
}  // namespace Geometry
}  // namespace Principia