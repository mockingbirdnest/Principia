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
  public T all { get; private set; }

  protected Boxed(T all) {
    this.all = all;
  }
}

// The |MarshalAsAttribute| does not support marshaling of generic types.
internal class BoxedInt32 : Boxed<int> {
  public static implicit operator BoxedInt32(int all) {
    return new BoxedInt32(all);
  }

  protected BoxedInt32(int all) : base(all) {}
}

}  // namespace ksp_plugin_adapter
}  // namespace principia
