using System;
using System.Collections.Generic;
using System.Reflection;
using System.Runtime.InteropServices;

namespace lamont {

// This class provides the following extension methods on all objects:
// — obj.Call("name")(args).
// — obj.GetValue("name");
// — obj.SetValue("name", value);
// The following generics are equivalent to casting the result of the
// non-generic versions, with better error messages:
// — obj.Call<T>("name")(args).
// — obj.GetValue<T>("name");
public static class Reflection {
  // Returns the value of the property or field of |obj| with the given name.
  public static T GetValue<T>(this object obj, string name) {
    if (obj == null) {
      throw new NullReferenceException(
          $"Cannot access {typeof(T).FullName} {name} on null object");
    }
    Type type = obj.GetType();
    object result = null;
    FieldInfo field = type.GetField(name, public_instance);
    PropertyInfo property = type.GetProperty(name, public_instance);
    if (field != null) {
      result = field.GetValue(obj);
    } else if (property != null) {
      result = property.GetValue(obj, index : null);
    } else {
      throw new MissingMemberException(
          $"No public instance field or property {name} in {type.FullName}");
    }
    try {
      return (T)result;
    } catch (Exception exception) {
      throw new InvalidCastException(
          $@"Could not convert the value of {
              (field == null ? "property" : "field")} {
              (field?.FieldType ?? property.PropertyType).FullName} {
              type.FullName}.{name}, {result}, to {typeof(T).FullName}",
          exception);
    }
  }

  public static void SetValue<T>(this object obj, string name, T value) {
    if (obj == null) {
      throw new NullReferenceException(
          $"Cannot set {typeof(T).FullName} {name} on null object");
    }
    Type type = obj.GetType();
    FieldInfo field = type.GetField(name, public_instance);
    PropertyInfo property = type.GetProperty(name, public_instance);
    if (field == null && property == null) {
      throw new MissingMemberException(
          $"No public instance field or property {name} in {type.FullName}");
    }
    try {
      field?.SetValue(obj, value);
      property?.SetValue(obj, value, index : null);
    } catch (Exception exception) {
      throw new ArgumentException(
          $@"Could not set {
              (field == null ? "property" : "field")} {
              (field?.FieldType ?? property.PropertyType).FullName} {
              type.FullName}.{name} to {typeof(T).FullName} {
              value?.GetType().FullName ?? "null"} {value}",
          exception);
    }
  }

  public static object GetValue(this object obj, string name) {
    return obj.GetValue<object>(name);
  }

  public delegate T BoundMethod<T>(params object[] args);

  public static BoundMethod<T> Call<T>(this object obj, string name) {
    if (obj == null) {
      throw new NullReferenceException($"Cannot call {name} on null object");
    }
    Type type = obj.GetType();
    MethodInfo method = type.GetMethod(name, public_instance);
    if (method == null) {
     throw new KeyNotFoundException(
         $"No public instance method {name} in {type.FullName}");
    }
    return args => {
      object result = method.Invoke(obj, args);
      try {
        return (T)result;
      } catch (Exception exception) {
        throw new InvalidCastException(
            $@"Could not convert the result of {
                method.ReturnType.FullName} {
                type.FullName}.{name}(), {result}, to {typeof(T).FullName}",
            exception);
      }
    };
  }

  public static BoundMethod<object> Call(this object obj, string name) {
    return obj.Call<object>(name);
  }

  private const BindingFlags public_instance =
      BindingFlags.Public | BindingFlags.Instance;
}

// Principia-specific utilities.
static class Principia {
  public static string AssemblyName() {
    foreach (var loaded_assembly in AssemblyLoader.loadedAssemblies) {
      if (loaded_assembly.assembly.GetName().Name == "principia.ksp_plugin_adapter") {
        return loaded_assembly.assembly.FullName;
      }
    }
    throw new DllNotFoundException(
        "principia.ksp_plugin_adapter not in AssemblyLoader.loadedAssemblies");
  }

  public static Type GetType(string name) {
    return Type.GetType(
      $"principia.ksp_plugin_adapter.{name}, {AssemblyName()}");
  }

  public static object Make(string type) {
    return Activator.CreateInstance(GetType(type));
  }

  public static object Get() {
    return GetType("ExternalInterface")
        .GetMethod("Get")
        .Invoke(null, null);
  }
}

[KSPAddon(startup : KSPAddon.Startup.EveryScene,
          once    : false)]
class LamontTest : ScenarioModule {
  private UnityEngine.Rect main_window_rectangle_;

  private void OnGUI() {
    main_window_rectangle_ = UnityEngine.GUILayout.Window(
        id         : this.GetHashCode(),
        screenRect : main_window_rectangle_,
        func       : DrawLamontWindow,
        text       : "LAMONT TEST",
        options    : UnityEngine.GUILayout.MinWidth(500));
  }

  void DrawLamontWindow(int window_id) {
    UnityEngine.GUILayout.BeginVertical();
    object principia = null;
    try {
      principia = Principia.Get();
    } catch (Exception e) {
      UnityEngine.GUILayout.TextArea(e.ToString());
    }
    try {
      var coefficient = principia.Call("GeopotentialGetCoefficient")(
          FlightGlobals.GetHomeBody().flightGlobalsIndex,
          2,
          0);
      var c20 = coefficient.GetValue<double>("x");
      UnityEngine.GUILayout.TextArea($"C20 = {c20}");
    } catch (Exception e) {
      UnityEngine.GUILayout.TextArea(e.ToString());
    }
    try {
      double r = principia.Call<double>("GeopotentialGetReferenceRadius")(
          FlightGlobals.GetHomeBody().flightGlobalsIndex);
      UnityEngine.GUILayout.TextArea($"r = {r}");
    } catch (Exception e) {
      UnityEngine.GUILayout.TextArea(e.ToString());
    }
    try {
      var xy = Principia.Make("XY");
      var c20 = xy.GetValue<string>("x");
    } catch (Exception e) {
      UnityEngine.GUILayout.TextArea(e.ToString());
    }
    try {
      var xy = Principia.Make("XY");
      var c20 = xy.GetValue<string>("z");
    } catch (Exception e) {
      UnityEngine.GUILayout.TextArea(e.ToString());
    }
    try {
      var xy = Principia.Make("Burn");
      var c20 = xy.GetValue<bool>("is_inertially_fixed");
    } catch (Exception e) {
      UnityEngine.GUILayout.TextArea(e.ToString());
    }
    try {
      var xy = Principia.Make("Burn");
      var c20 = xy.GetValue<string>("is_inertially_fixed");
    } catch (Exception e) {
      UnityEngine.GUILayout.TextArea(e.ToString());
    }
    try {
      var qp = Principia.Make("QP");
      var q = qp.GetValue("q");
      var x = q.GetValue<double>("x");
      q.SetValue("y", "kitten");
    } catch (Exception e) {
      UnityEngine.GUILayout.TextArea(e.ToString());
    }
    UnityEngine.GUILayout.EndVertical();
  }
}

}
