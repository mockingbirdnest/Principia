using System;
using System.Collections.Generic;
using System.Runtime.InteropServices;

namespace principia {
namespace ksp_plugin_adapter {

[KSPAddon(KSPAddon.Startup.Flight, false)]
public class PluginAdapter : UnityEngine.MonoBehaviour {

  #region Interface
  private const string kDllPath = "GameData/Principia/principia.dll";

  [StructLayout(LayoutKind.Sequential)]
  private struct XYZ {
    double x, y, z;
    public static explicit operator XYZ(Vector3d v) {
      return new XYZ {x = v.x, y = v.y, z = v.z};
    }
    public static explicit operator Vector3d(XYZ v) {
      return new Vector3d {x = v.x, y = v.y, z = v.z};
    }
  };

  [DllImport(dllName           : kDllPath,
             CallingConvention = CallingConvention.Cdecl)]
  private static extern void InitGoogleLogging();

  [DllImport(dllName           : kDllPath,
             CallingConvention = CallingConvention.Cdecl)]
  private static extern IntPtr SayHello();

  [DllImport(dllName           : kDllPath,
             CallingConvention = CallingConvention.Cdecl)]
  private static extern IntPtr CreatePlugin(
      double initial_time,
      int sun_index,
      double sun_gravitational_parameter,
      double planetarium_rotation_in_degrees);

  [DllImport(dllName           : kDllPath,
             CallingConvention = CallingConvention.Cdecl)]
  private static extern void DestroyPlugin(IntPtr plugin);

  [DllImport(dllName           : kDllPath,
             CallingConvention = CallingConvention.Cdecl)]
  private static extern void InsertCelestial(
      IntPtr plugin,
      int index,
      double gravitational_parameter,
      int parent,
      XYZ from_parent_position,
      XYZ from_parent_velocity);

  [DllImport(dllName           : kDllPath,
             CallingConvention = CallingConvention.Cdecl)]
  private static extern void UpdateCelestialHierarchy(IntPtr plugin, int index,
                                                      int parent);

  [DllImport(dllName           : kDllPath,
             CallingConvention = CallingConvention.Cdecl)]
  private static extern bool InsertOrKeepVessel(
      IntPtr plugin,
      [MarshalAs(UnmanagedType.LPStr)]
      String guid,
      int parent);

  [DllImport(dllName           : kDllPath,
             CallingConvention = CallingConvention.Cdecl)]
  private static extern void SetVesselStateOffset(
      IntPtr plugin,
      [MarshalAs(UnmanagedType.LPStr)]
      String guid,
      XYZ from_parent_position,
      XYZ from_parent_velocity);

  [DllImport(dllName           : kDllPath,
             CallingConvention = CallingConvention.Cdecl)]
  private static extern XYZ VesselDisplacementFromParent(
      IntPtr plugin,
      [MarshalAs(UnmanagedType.LPStr)]
      String guid);

  [DllImport(dllName           : kDllPath,
             CallingConvention = CallingConvention.Cdecl)]
  private static extern XYZ VesselParentRelativeVelocity(
      IntPtr plugin,
      [MarshalAs(UnmanagedType.LPStr)]
      String guid);

  [DllImport(dllName           : kDllPath,
             CallingConvention = CallingConvention.Cdecl)]
  private static extern XYZ CelestialDisplacementFromParent(
      IntPtr plugin,
      int index);
  
  [DllImport(dllName           : kDllPath,
             CallingConvention = CallingConvention.Cdecl)]
  private static extern XYZ CelestialParentRelativeVelocity(
      IntPtr plugin,
      int index);
  #endregion

  private UnityEngine.Rect window_position_;
  private IntPtr plugin_;
  private bool plugin_running_;

  ~PluginAdapter() {
    DestroyPlugin(plugin_);
  }

  private delegate void BodyProcessor(CelestialBody body);
  private delegate void VesselProcessor(Vessel vessel);

  // Applies |process_body| to all bodies but the sun in the tree of celestial
  // bodies, in topological order.
  private void ApplyToBodyTree(BodyProcessor process_body) {
    // Tree traversal (DFS, not that it matters).
    Stack<CelestialBody> stack = new Stack<CelestialBody>();
    foreach (CelestialBody child in Planetarium.fetch.Sun.orbitingBodies) {
      stack.Push(child);
    }
    CelestialBody body;
    while (stack.Count > 0) {
      body = stack.Pop();
      process_body(body);
      foreach (CelestialBody child in body.orbitingBodies) {
        stack.Push(child);
      }
    }
  }

  // Applies |process_vessel| to all vessels in space.
  private void ApplyToVesselsInSpace(VesselProcessor process_vessel) {
    foreach (Vessel vessel in FlightGlobals.Vessels) {
      if (vessel.situation == Vessel.Situations.SUB_ORBITAL ||
          vessel.situation == Vessel.Situations.ORBITING ||
          vessel.situation == Vessel.Situations.ESCAPING) {
        process_vessel(vessel);
      }
    }
  }

  #region Unity Lifecycle
  private void Start() {
    InitGoogleLogging();
    RenderingManager.AddToPostDrawQueue(queueSpot    : 3,
                                        drawFunction : new Callback(DrawGUI));
    window_position_ = new UnityEngine.Rect(
        left   : UnityEngine.Screen.width / 2.0f,
        top    : UnityEngine.Screen.height / 2.0f,
        width  : 10,
        height : 10);
  }

  private void OnDestroy() {
    RenderingManager.RemoveFromPostDrawQueue(
        queueSpot    : 3,
        drawFunction : new Callback(DrawGUI));
  }

  private void FixedUpdate() {
    if (plugin_running_) {
      double universal_time = Planetarium.GetUniversalTime();
      BodyProcessor update_body = body => {
        Vector3d position =
            (Vector3d)CelestialDisplacementFromParent(plugin_,
                                                      body.flightGlobalsIndex);
        Vector3d velocity =
            (Vector3d)CelestialParentRelativeVelocity(plugin_,
                                                      body.flightGlobalsIndex);
        body.orbit.UpdateFromStateVectors(pos: position, vel: velocity,
                                          refBody: body.orbit.referenceBody,
                                          UT: universal_time);
      };
      ApplyToBodyTree(update_body);
      VesselProcessor update_vessel = vessel => {
        Vector3d position =
            (Vector3d)VesselDisplacementFromParent(plugin_,
                                                   vessel.id.ToString());
        Vector3d velocity =
            (Vector3d)VesselParentRelativeVelocity(plugin_,
                                                   vessel.id.ToString());
        vessel.orbit.UpdateFromStateVectors(pos: position, vel: velocity,
                                            refBody: vessel.orbit.referenceBody,
                                            UT: universal_time);
      };
      ApplyToVesselsInSpace(update_vessel);
    }
  }
  #endregion

  private void DrawGUI() {
    UnityEngine.GUI.skin = HighLogic.Skin;
    window_position_ = UnityEngine.GUILayout.Window(
        id         : 1,
        screenRect : window_position_,
        func       : DrawMainWindow,
        text       : "Traces of Various Descriptions",
        options    : UnityEngine.GUILayout.MinWidth(500));
  }

  private void DrawMainWindow(int window_id) {
    UnityEngine.GUIStyle style = new UnityEngine.GUIStyle(
        UnityEngine.GUI.skin.button);
    style.normal.textColor = style.focused.textColor = UnityEngine.Color.white;
    style.hover.textColor = style.active.textColor = UnityEngine.Color.yellow;
    style.onNormal.textColor  = UnityEngine.Color.green;
    style.onFocused.textColor = UnityEngine.Color.green;
    style.onHover.textColor   = UnityEngine.Color.green;
    style.onActive.textColor  = UnityEngine.Color.green;
    style.padding             = new UnityEngine.RectOffset(8, 8, 8, 8);

    UnityEngine.GUILayout.BeginVertical();
    IntPtr hello_ptr = SayHello();
    UnityEngine.GUILayout.TextArea(text : Marshal.PtrToStringAnsi(hello_ptr));
    if (UnityEngine.GUILayout.Button(plugin_running_? "Stop plugin"
                                                    : "Start plugin")) {
      if (plugin_running_) {
        DestroyPlugin(plugin_);
      } else {
        InitializePlugin();
      }
      plugin_running_ = !plugin_running_;
    }
    UnityEngine.GUILayout.EndVertical();

    UnityEngine.GUI.DragWindow(
        position : new UnityEngine.Rect(left : 0f, top : 0f, width : 10000f,
                                        height : 20f));
  }

  private void InitializePlugin() {;
    plugin_ = CreatePlugin(Planetarium.GetUniversalTime(),
                           Planetarium.fetch.Sun.flightGlobalsIndex,
                           Planetarium.fetch.Sun.gravParameter,
                           Planetarium.InverseRotAngle);
    ApplyToBodyTree(
        body => InsertCelestial(plugin_, body.flightGlobalsIndex,
                                body.gravParameter,
                                body.orbit.referenceBody.flightGlobalsIndex,
                                (XYZ)body.orbit.pos, (XYZ)body.orbit.vel));
    VesselProcessor insert_vessel = vessel => {
      bool inserted =
          InsertOrKeepVessel(plugin_, vessel.id.ToString(),
                             vessel.orbit.referenceBody.flightGlobalsIndex);
      if (!inserted) {
        UnityEngine.Debug.LogError(
            "Plugin initialisation: vessel not inserted");
      }
    };
    ApplyToVesselsInSpace(insert_vessel);
  }
}

}  // namespace ksp_plugin_adapter
}  // namespace principia
