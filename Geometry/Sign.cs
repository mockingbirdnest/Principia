using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Geometry {
  public struct Sign {
    public readonly bool Positive;
    private Sign(bool positive) { Positive = positive; }
    public static explicit operator Sign(int x) {
      return new Sign(x > 0);
    }
    public static explicit operator Sign(Scalar x) {
      return new Sign(x > (Scalar)0);
    }
    public static implicit operator Scalar(Sign sign) {
      return sign.Positive ? (Scalar)1 : (Scalar)(-1);
    }
    public static Sign operator *(Sign left, Sign right) {
      return new Sign(left.Positive == right.Positive);
    }
  }
}