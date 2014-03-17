#include "stdafx.h"

using namespace System;
using namespace System::Text;
using namespace System::Collections::Generic;
using namespace Microsoft::VisualStudio::TestTools::UnitTesting;
using namespace Geometry;

namespace GeometryTests {
  [TestClass]
  public ref class UnitTest {
  private:
    TestContext^ testContextInstance;

    ref class A : ISpace {};

#pragma warning(disable:4395)
    double const ε = 1e-13;
    template<typename Vector>
    void AssertNullVector(Vector v) {
      Assert::AreEqual((double)v.Coordinates.X, 0.0, ε);
      Assert::AreEqual((double)v.Coordinates.Y, 0.0, ε);
      Assert::AreEqual((double)v.Coordinates.Z, 0.0, ε);
    }
    template<typename Vector>
    void AssertNotNullVector(Vector v) {
      Assert::AreNotEqual((double)v.Coordinates.X, 0.0, ε);
      Assert::AreNotEqual((double)v.Coordinates.Y, 0.0, ε);
      Assert::AreNotEqual((double)v.Coordinates.Z, 0.0, ε);
    }
#pragma warning(default:4395)

    template<typename Vector>
    void AssertVectorsEqual(Vector v, Vector w) {
      AssertNullVector(v - w);
    }

    template<typename Vector>
    void AssertVectorsNotEqual(Vector v, Vector w) {
      AssertNotNullVector(v - w);
    }

    template<typename Vector>
    void TestInnerProductSpaceOperations() {
      Vector u = Vector(R3Element((Scalar)9.1, (Scalar)(-0.1), (Scalar)10.0));
      Vector t = Vector(R3Element((Scalar)0.0, (Scalar)1.9, (Scalar)(-0.5)));
      Vector v = Vector(R3Element((Scalar)1.0, (Scalar)(-2.0), (Scalar)3.0));
      Vector w = Vector(R3Element((Scalar)4.0, (Scalar)5.0, (Scalar)0.5));
      Vector i = Vector(R3Element((Scalar)1.0, (Scalar)0.0, (Scalar)0.0));
      Vector j = Vector(R3Element((Scalar)0.0, (Scalar)1.0, (Scalar)0.0));
      Vector k = Vector(R3Element((Scalar)0.0, (Scalar)0.0, (Scalar)1.0));
      Scalar α = (Scalar)42.0;
      Scalar β = (Scalar)(-6.0 * 9.0);
      // Nontriviality.
      AssertNotNullVector(v);
      // Commutativity of addition.
      AssertVectorsEqual(v + w, w + v);
      // Additive inverse and subtraction.
      AssertNullVector(v - v);
      AssertVectorsEqual(v - w, v + (-w));
      AssertVectorsEqual(v - w, -(w - v));
      // Noncommutativity of subtraction.
      AssertVectorsNotEqual(v - w, w - v);
      // Multiplication by a scalar.
      AssertVectorsEqual(v + v, (Scalar)2.0 * v);
      AssertVectorsEqual(α * v, v * α);
      // Multiplicative inverse of scalars.
      AssertVectorsEqual((α * v) / α, v);
      // Bilinearity of the inner product.
      Assert::AreEqual(Vector::InnerProduct(α * u + v, β * w + t),
                       α * β * Vector::InnerProduct(u, w)
                       + α * Vector::InnerProduct(u, t)
                       + β * Vector::InnerProduct(v, w)
                       + Vector::InnerProduct(v, t));
      // Commutativity of the inner product.
      Assert::AreEqual(Vector::InnerProduct(v, w), Vector::InnerProduct(w, v));
      // Inner product scaling.
      Assert::AreEqual(Vector::InnerProduct(j, j), (Scalar)1.0);
      // Standard basis orthogonality.
      Assert::AreEqual(Vector::InnerProduct(i, j), (Scalar)0.0);
      Assert::AreEqual(Vector::InnerProduct(j, k), (Scalar)0.0);
      Assert::AreEqual(Vector::InnerProduct(k, i), (Scalar)0.0);
    }

    void TestSO3() {
      Vector<A^> v = Vector<A^>(R3Element((Scalar)1.0, (Scalar)(-2.0), (Scalar)3.0));
      BiVector<A^> x = BiVector<A^>(R3Element((Scalar)1.0, (Scalar)(-3.0), (Scalar) 2.0));
      BiVector<A^> y = BiVector<A^>(R3Element((Scalar)6.0, (Scalar)(-2.7), (Scalar)5.0));
      BiVector<A^> z = BiVector<A^>(R3Element((Scalar)(-10.0), (Scalar)0.0, (Scalar)(-.025)));
      Scalar λ = (Scalar)12.345;
      // Nontriviality of the commutator.
      AssertNotNullVector(BiVector<A^>::Commutator(x, y));
      // Linearity of the commutator in the first argument.
      AssertVectorsEqual(BiVector<A^>::Commutator(λ * x + z, y),
                         λ * BiVector<A^>::Commutator(x, y)
                         + BiVector<A^>::Commutator(z, y));
      // Anticommutativity of the commutator.
      AssertVectorsEqual(BiVector<A^>::Commutator(x, y),
                         -BiVector<A^>::Commutator(y, x));
      AssertNullVector(BiVector<A^>::Commutator(x, x));
      // Jacobi identity.
      AssertNullVector(
        BiVector<A^>::Commutator(x, BiVector<A^>::Commutator(y, z))
        + BiVector<A^>::Commutator(y, BiVector<A^>::Commutator(z, x))
        + BiVector<A^>::Commutator(z, BiVector<A^>::Commutator(x, y)));
      // Left and right actions.
      AssertVectorsEqual(v.ActedUponBy(x), -x.ActOn(v));
      AssertVectorsNotEqual(v.ActedUponBy(x), x.ActOn(v));
      // Orthogonality of the commutator.
      Assert::AreEqual(
        BiVector<A^>::InnerProduct(BiVector<A^>::Commutator(x, y), x),
        (Scalar)0.0);
      Assert::AreEqual(
        BiVector<A^>::InnerProduct(BiVector<A^>::Commutator(x, y), y),
        (Scalar)0.0);
      // Exponential map.
      AssertVectorsEqual(
        BiVector<A^>::Exp(
          BiVector<A^>(R3Element((Scalar)0.0, (Scalar)0.0, (Scalar)0.0))
          ).ActOn(v),
        v);
      AssertVectorsEqual(BiVector<A^>::Exp(-x).ActOn(v),
                         BiVector<A^>::Exp(x).Inverse().ActOn(v));
      AssertVectorsEqual(
        BiVector<A^>::Exp(x + λ * x).ActOn(v),
        Maps::Compose(BiVector<A^>::Exp(x), BiVector<A^>::Exp(λ * x)).ActOn(v));
    }

  public:
    /// <summary>
    ///Gets or sets the test context which provides
    ///information about and functionality for the current test run.
    ///</summary>
    property Microsoft::VisualStudio::TestTools::UnitTesting::TestContext^ TestContext
    {
      Microsoft::VisualStudio::TestTools::UnitTesting::TestContext^ get() {
        return testContextInstance;
      }
      System::Void set(Microsoft::VisualStudio::TestTools::UnitTesting::TestContext^ value) {
        testContextInstance = value;
      }
    };

#pragma region Additional test attributes
    //
    //You can use the following additional attributes as you write your tests:
    //
    //Use ClassInitialize to run code before running the first test in the class
    //[ClassInitialize()]
    //static void MyClassInitialize(TestContext^ testContext) {};
    //
    //Use ClassCleanup to run code after all tests in a class have run
    //[ClassCleanup()]
    //static void MyClassCleanup() {};
    //
    //Use TestInitialize to run code before running each test
    //[TestInitialize()]
    //void MyTestInitialize() {};
    //
    //Use TestCleanup to run code after each test has run
    //[TestCleanup()]
    //void MyTestCleanup() {};
    //
#pragma endregion

    [TestMethod]
    void GrassmannAlgebra() {
      TestInnerProductSpaceOperations<Vector<A^>>();
      TestInnerProductSpaceOperations<BiVector<A^>>();
      TestSO3();
    };
  };
}