using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Geometry {
  public struct Sign {
    public bool positive;
    public static explicit operator Sign(int x) {
      return new Sign { positive = x > 0 };
    }
    public static explicit operator Sign(Scalar x) {
      return new Sign { positive = x > (Scalar)0 };
    }
    public static implicit operator Scalar(Sign sgn) {
      return sgn.positive ? (Scalar)1 : (Scalar)(-1);
    }
    public static Sign operator *(Sign left, Sign right) {
      return new Sign { positive = (left.positive == right.positive) };
    }
  }
}