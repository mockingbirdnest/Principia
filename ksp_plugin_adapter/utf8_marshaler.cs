using System;
using System.Globalization;
using System.Runtime.InteropServices;
using System.Text;

namespace principia {
namespace ksp_plugin_adapter {

internal abstract class UTF8Marshaler : ICustomMarshaler {
  public abstract void CleanUpNativeData(IntPtr native_data);
  public abstract IntPtr MarshalManagedToNative(object managed_object);
  public abstract object MarshalNativeToManaged(IntPtr native_data);

  void ICustomMarshaler.CleanUpManagedData(object managed_object) {}

  int ICustomMarshaler.GetNativeDataSize() {
    return -1;
  }

  protected readonly static Encoding utf8_ =
      new UTF8Encoding(encoderShouldEmitUTF8Identifier : false,
                       throwOnInvalidBytes             : true);
}

// A marshaler for in parameter UTF-8 strings whose ownership is not taken from
// the caller.
internal class InUTF8Marshaler : UTF8Marshaler {
  // In addition to implementing the |ICustomMarshaler| interface, custom
  // marshalers must implement a static method called |GetInstance| that accepts
  // a |String| as a parameter and has a return type of |ICustomMarshaler|,
  // see https://goo.gl/wwmBTa.
  public static ICustomMarshaler GetInstance(String s) {
    return instance_;
  }

  public override void CleanUpNativeData(IntPtr native_data) {
    Marshal.FreeHGlobal(native_data);
#if MARSHAL_DEBUG
    // Useful for debugging things from a console application to see whether
    // allocating and freeing match.  We cannot use Log.Info here, since that
    // would recursively call the marshaler and end in overflowing hilarity.
    Console.WriteLine("freeh " + Convert.ToString(native_data.ToInt64(), 16));
#endif
  }

  public override IntPtr MarshalManagedToNative(object managed_object) {
    var value = managed_object as String;
    if (value == null) {
      throw Log.Fatal(String.Format(CultureInfo.InvariantCulture,
                                    "|{0}| must be used on a |{1}|.",
                                    GetType().Name,
                                    typeof(String).Name));
    }
    int size = utf8_.GetByteCount(value);
    IntPtr ptr = Marshal.AllocHGlobal(size + 1);
    unsafe {
      byte* begin = (byte*)ptr;
#if MARSHAL_DEBUG
      Console.WriteLine("alloch " + Convert.ToString(ptr.ToInt64(), 16));
#endif
      fixed (char* utf16_units = value) {
        utf8_.GetBytes(utf16_units, value.Length, begin, size);
      }
      begin[size] = 0;
    }
    return ptr;
  }

  public override object MarshalNativeToManaged(IntPtr native_data) {
    throw Log.Fatal("use |OutUTF8Marshaler| for out parameters");
  }

  private readonly static InUTF8Marshaler instance_ = new InUTF8Marshaler();
}

// A marshaler for out parameter or return value UTF-8 strings whose ownership
// is not taken by the caller.
internal class OutUTF8Marshaler : UTF8Marshaler {
  public static ICustomMarshaler GetInstance(String s) {
    return instance_;
  }

  public override void CleanUpNativeData(IntPtr native_data) {}
  public override IntPtr MarshalManagedToNative(object managed_object) {
    throw Log.Fatal("use |InUTF8Marshaler| for in parameters");
  }

  public override unsafe object MarshalNativeToManaged(IntPtr native_data) {
    sbyte* begin = (sbyte*)native_data;
    sbyte* end;
    for (end = begin; *end != 0; ++end) {}
    return new String(begin, 0, (int)(end - begin), utf8_);
  }

  private readonly static OutUTF8Marshaler instance_ = new OutUTF8Marshaler();
}

// A marshaler for out parameter or return value UTF-8 strings whose ownership
// is taken by the caller.
internal class UTF8FactoryMarshaler : OutUTF8Marshaler {
  public static new ICustomMarshaler GetInstance(String s) {
    return instance_;
  }

  public override void CleanUpNativeData(IntPtr native_data) {
    // NOTE(egg): this should get renamed.
    Interface.DeletePluginSerialization(ref native_data);
  }

  private readonly static UTF8FactoryMarshaler instance_ =
      new UTF8FactoryMarshaler();
}

}  // namespace ksp_plugin_adapter
}  // namespace principia
