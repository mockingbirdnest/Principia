using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;


namespace principia {
namespace ksp_plugin_adapter {

// Strongly-typed boxing.  This is used for marshaling, since custom marshalers
// are only allowed on classes, strings, arrays, and boxed value types, so that
// |T?| cannot be marshaled, and |object| is not statically typed.
internal class Box<T> where T : struct {
  public static implicit operator Box<T>(T all) {
    return new Box<T>(all);
  }

  public T all { get; }

  private Box(T all) {
    this.all = all;
  }
}

}  // namespace ksp_plugin_adapter
}  // namespace principia
