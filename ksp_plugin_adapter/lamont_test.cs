using System;
using System.Reflection;

namespace lamont {

static class Reflection {
  public class ReflectedProperties {
    public ReflectedProperties(object obj) {
      obj_ = obj;
    }

    public object this[string name] {
      get {
        return obj_.GetType()
            .GetProperty(name, BindingFlags.Public | BindingFlags.Instance)
            .GetValue(obj_, index : null);
      }
      set {
        obj_.GetType()
            .GetProperty(name, BindingFlags.Public | BindingFlags.Instance)
            .SetValue(obj_, value, index : null);
      }
    }

    private object obj_;
  }

  public class ReflectedFields {
    public ReflectedFields(object obj) {
      obj_ = obj;
    }

    public object this[string name] {
      get {
        return obj_.GetType()
            .GetField(name, BindingFlags.Public | BindingFlags.Instance)
            .GetValue(obj_);
      }
      set {
        obj_.GetType()
            .GetField(name, BindingFlags.Public | BindingFlags.Instance)
            .SetValue(obj_, value);
      }
    }

    private object obj_;
  }

  public class ReflectedNonPublicFields {
    public ReflectedNonPublicFields(object obj) {
      obj_ = obj;
    }

    public object this[string name] {
      get {
        return obj_.GetType()
            .GetField(name, BindingFlags.NonPublic | BindingFlags.Instance)
            .GetValue(obj_);
      }
      set {
        obj_.GetType()
            .GetField(name, BindingFlags.NonPublic | BindingFlags.Instance)
            .SetValue(obj_, value);
      }
    }

    private object obj_;
  }

  public static ReflectedProperties Properties(this object obj){
     return new ReflectedProperties(obj);
  }

  public static ReflectedFields Fields(this object obj){
     return new ReflectedFields(obj);
  }

  public static ReflectedNonPublicFields NonPublicFields(this object obj){
     return new ReflectedNonPublicFields(obj);
  }
}

static class Principia {
  public class Error : Exception {
    public Error() {}
    public Error(string message) : base(message) {}
    Error(string message, Exception inner) : base(message, inner) {}
  }

  public static string AssemblyName() {
    foreach (var loaded_assembly in AssemblyLoader.loadedAssemblies) {
      if (loaded_assembly.assembly.GetName().Name == "ksp_plugin_adapter") {
        return loaded_assembly.assembly.FullName;
      }
    }
    throw new DllNotFoundException(
        "ksp_plugin_adapter not in AssemblyLoader.loadedAssemblies");
  }

  public static Type GetType(string name) {
    return Type.GetType(
      $"principia.ksp_plugin_adapter.{name}, {AssemblyName()}");
  }

  public static object Make(string name) {
    return Activator.CreateInstance(GetType(name));
  }

  public delegate void InterfaceFunction(params object[] args);

  public static InterfaceFunction Call(string name) {
    return args => {
      object status = GetType("Interface").GetMethod(
          $"External{name}",
          BindingFlags.NonPublic | BindingFlags.Static).Invoke(null, args);
      if ((int)status.Fields()["error"] != 0) {
        throw new Error(
          $"Error {status.Fields()["error"]} from call to {name} with arguments {args}");
      }
    };
  }

  public static IntPtr Plugin() {
    foreach (var module in ScenarioRunner.GetLoadedModules()) {
      if (module.GetType().FullName ==
          "principia.ksp_plugin_adapter.PrincipiaPluginAdapter") {
        return (IntPtr)module.NonPublicFields()["plugin_"];
      }
    }
    return IntPtr.Zero;
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
    try {
      object coefficient = Principia.Make("XY");
      Principia.Call("GeopotentialGetCoefficient")(
          Principia.Plugin(),
          FlightGlobals.GetHomeBody().flightGlobalsIndex,
          2,
          0,
          coefficient);
      double c20 = (double)coefficient.Fields()["x"];
      UnityEngine.GUILayout.TextArea($"C20 = {c20}");
    } catch (Exception e) {
      UnityEngine.GUILayout.TextArea(e.ToString());
    }
    try {
      double r = double.NaN;
      Principia.Call("GeopotentialGetReferenceRadius")(
          Principia.Plugin(),
          FlightGlobals.GetHomeBody().flightGlobalsIndex,
          r);
      UnityEngine.GUILayout.TextArea($"r = {r}");
    } catch (Exception e) {
      UnityEngine.GUILayout.TextArea(e.ToString());
    }
    UnityEngine.GUILayout.EndVertical();
  }
}

}
