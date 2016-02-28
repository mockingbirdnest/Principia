using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;


namespace principia {
namespace ksp_plugin_adapter {

// Strongly-typed boxing.  This is used for marshaling optional parameters,
// since custom marshalers are only allowed on classes, strings, arrays, and
// boxed value types, so that |T?| cannot be marshaled, and |object| is not
// statically typed.
internal class Boxed<T> where T : struct {
  public static implicit operator Boxed<T>(T all) {
    return new Boxed<T>(all);
  }

  public T all { get; private set; }

  protected Boxed(T all) {
    this.all = all;
  }
}

internal class BoxedInteger : Boxed<int> {
  protected BoxedInteger(int all) : base(all) {}
}

}  // namespace ksp_plugin_adapter
}  // namespace principia
