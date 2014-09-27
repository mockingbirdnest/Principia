using System;
using System.Collections.Generic;
using System.Runtime.InteropServices;

namespace principia {
namespace ksp_plugin_adapter {

[KSPAddon(startup : KSPAddon.Startup.Flight, once : false)]
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

  #region glog
  [DllImport(dllName           : kDllPath,
             CallingConvention = CallingConvention.Cdecl)]
  private static extern void InitGoogleLogging();

  [DllImport(dllName           : kDllPath,
             CallingConvention = CallingConvention.Cdecl)]
  private static extern void LogInfo(
      [MarshalAs(UnmanagedType.LPStr)] String message);

  [DllImport(dllName           : kDllPath,
             CallingConvention = CallingConvention.Cdecl)]
  private static extern void LogWarning(
      [MarshalAs(UnmanagedType.LPStr)] String message);

  [DllImport(dllName           : kDllPath,
             CallingConvention = CallingConvention.Cdecl)]
  private static extern void LogError(
      [MarshalAs(UnmanagedType.LPStr)] String message);

  [DllImport(dllName           : kDllPath,
             CallingConvention = CallingConvention.Cdecl)]
  private static extern void LogFatal(
      [MarshalAs(UnmanagedType.LPStr)] String message);
  #endregion

  [DllImport(dllName           : kDllPath,
             CallingConvention = CallingConvention.Cdecl)]
  private static extern IntPtr SayHello();

  [DllImport(dllName           : kDllPath,
             CallingConvention = CallingConvention.Cdecl)]
  private static extern IntPtr NewPlugin(
      double initial_time,
      int sun_index,
      double sun_gravitational_parameter,
      double planetarium_rotation_in_degrees);

  [DllImport(dllName           : kDllPath,
             CallingConvention = CallingConvention.Cdecl)]
  private static extern void DeletePlugin(ref IntPtr plugin);

  [DllImport(dllName           : kDllPath,
             CallingConvention = CallingConvention.Cdecl)]
  private static extern void InsertCelestial(
      IntPtr plugin,
      int celestial_index,
      double gravitational_parameter,
      int parent_index,
      XYZ from_parent_position,
      XYZ from_parent_velocity);

  [DllImport(dllName           : kDllPath,
             CallingConvention = CallingConvention.Cdecl)]
  private static extern void UpdateCelestialHierarchy(IntPtr plugin,
                                                      int celestial_index,
                                                      int parent_index);

  [DllImport(dllName           : kDllPath,
             CallingConvention = CallingConvention.Cdecl)]
  private static extern bool InsertOrKeepVessel(
      IntPtr plugin,
      [MarshalAs(UnmanagedType.LPStr)] String guid,
      int parent_index);

  [DllImport(dllName           : kDllPath,
             CallingConvention = CallingConvention.Cdecl)]
  private static extern void SetVesselStateOffset(
      IntPtr plugin,
      [MarshalAs(UnmanagedType.LPStr)] String guid,
      XYZ from_parent_position,
      XYZ from_parent_velocity);

  [DllImport(dllName           : kDllPath,
             CallingConvention = CallingConvention.Cdecl)]
  private static extern void AdvanceTime(IntPtr plugin, 
                                         double t,
                                         double planetarium_rotation);

  [DllImport(dllName           : kDllPath,
             CallingConvention = CallingConvention.Cdecl)]
  private static extern XYZ VesselDisplacementFromParent(
      IntPtr plugin,
      [MarshalAs(UnmanagedType.LPStr)] String guid);

  [DllImport(dllName           : kDllPath,
             CallingConvention = CallingConvention.Cdecl)]
  private static extern XYZ VesselParentRelativeVelocity(
      IntPtr plugin,
      [MarshalAs(UnmanagedType.LPStr)] String guid);

  [DllImport(dllName           : kDllPath,
             CallingConvention = CallingConvention.Cdecl)]
  private static extern XYZ CelestialDisplacementFromParent(
      IntPtr plugin,
      int celestial_index);
  
  [DllImport(dllName           : kDllPath,
             CallingConvention = CallingConvention.Cdecl)]
  private static extern XYZ CelestialParentRelativeVelocity(
      IntPtr plugin,
      int celestial_index);
  #endregion

  private UnityEngine.Rect window_position_;
  private IntPtr plugin_ = IntPtr.Zero;

  PluginAdapter() {
    // We create this directory here so we do not need to worry about cross-
    // platform problems in C++.
    System.IO.Directory.CreateDirectory("glog/Principia");
    InitGoogleLogging();
  }

  ~PluginAdapter() {
    DeletePlugin(ref plugin_);
  }

  private bool PluginRunning() {
    return plugin_ != IntPtr.Zero;
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
  // See the Unity manual on execution order for more information on |Start()|,
  // |OnDestroy()| and |FixedUpdate()|.
  // http://docs.unity3d.com/Manual/ExecutionOrder.html
  private void Start() {
    LogInfo("principia.ksp_plugin_adapter.PluginAdapter.Start()");
    RenderingManager.AddToPostDrawQueue(queueSpot    : 3,
                                        drawFunction : new Callback(DrawGUI));
    window_position_ = new UnityEngine.Rect(
        left   : UnityEngine.Screen.width / 2.0f,
        top    : UnityEngine.Screen.height / 2.0f,
        width  : 10,
        height : 10);
  }

  private void OnDestroy() {
    LogInfo("principia.ksp_plugin_adapter.PluginAdapter.OnDestroy()");
    RenderingManager.RemoveFromPostDrawQueue(
        queueSpot    : 3,
        drawFunction : new Callback(DrawGUI));
  }

  private void FixedUpdate() {
    if (PluginRunning()) {
      double universal_time = Planetarium.GetUniversalTime();
      AdvanceTime(plugin_, universal_time, Planetarium.InverseRotAngle);
      BodyProcessor update_body = body => {
        UpdateCelestialHierarchy(plugin_,
                                 body.flightGlobalsIndex,
                                 body.orbit.referenceBody.flightGlobalsIndex);
        Vector3d position =
            (Vector3d)CelestialDisplacementFromParent(plugin_,
                                                      body.flightGlobalsIndex);
        Vector3d velocity =
            (Vector3d)CelestialParentRelativeVelocity(plugin_,
                                                      body.flightGlobalsIndex);
        // TODO(egg): Some of this might be be superfluous and redundant.
        Orbit original = body.orbit;
        Orbit copy = new Orbit(original.inclination, original.eccentricity,
                               original.semiMajorAxis, original.LAN,
                               original.argumentOfPeriapsis,
                               original.meanAnomalyAtEpoch, original.epoch,
                               original.referenceBody);
        copy.UpdateFromStateVectors(position, velocity, copy.referenceBody,
                                    universal_time);
        body.orbit.inclination = copy.inclination;
        body.orbit.eccentricity = copy.eccentricity;
        body.orbit.semiMajorAxis = copy.semiMajorAxis;
        body.orbit.LAN = copy.LAN;
        body.orbit.argumentOfPeriapsis = copy.argumentOfPeriapsis;
        body.orbit.meanAnomalyAtEpoch = copy.meanAnomalyAtEpoch;
        body.orbit.epoch = copy.epoch;
        body.orbit.referenceBody = copy.referenceBody;;
        body.orbit.Init();
        body.orbit.UpdateFromUT(universal_time);
        body.CBUpdate();
        body.orbit.UpdateFromStateVectors((Vector3d)position,
                                          (Vector3d)velocity,
                                          copy.referenceBody, universal_time);
      };
      ApplyToBodyTree(update_body);
      VesselProcessor update_vessel = vessel => {
        InsertOrKeepVessel(plugin_,
                           vessel.id.ToString(),
                           vessel.orbit.referenceBody.flightGlobalsIndex);
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
    if (UnityEngine.GUILayout.Button(PluginRunning() ? "Stop plugin"
                                                     : "Start plugin")) {
      if (PluginRunning()) {
        DeletePlugin(ref plugin_);
      } else {
        InitializePlugin();
      }
    }
    UnityEngine.GUILayout.EndVertical();

    UnityEngine.GUI.DragWindow(
        position : new UnityEngine.Rect(left : 0f, top : 0f, width : 10000f,
                                        height : 20f));
  }

  private void InitializePlugin() {;
    plugin_ = NewPlugin(Planetarium.GetUniversalTime(),
                        Planetarium.fetch.Sun.flightGlobalsIndex,
                        Planetarium.fetch.Sun.gravParameter,
                        Planetarium.InverseRotAngle);
    BodyProcessor insert_body = body => {
      LogInfo("Inserting " + body.name + "...");
      InsertCelestial(plugin_, body.flightGlobalsIndex, body.gravParameter,
                      body.orbit.referenceBody.flightGlobalsIndex,
                      (XYZ)body.orbit.pos, (XYZ)body.orbit.vel);
    };
    ApplyToBodyTree(insert_body);
    VesselProcessor insert_vessel = vessel => {
      LogInfo("Inserting " + vessel.name + "...");
      bool inserted =
          InsertOrKeepVessel(plugin_, vessel.id.ToString(),
                             vessel.orbit.referenceBody.flightGlobalsIndex);
      if (!inserted) {
        LogFatal("Plugin initialisation: vessel not inserted");
      } else {
        SetVesselStateOffset(plugin_, vessel.id.ToString(),
                             (XYZ)vessel.orbit.pos, (XYZ)vessel.orbit.vel);
      }
    };
    ApplyToVesselsInSpace(insert_vessel);
  }
}

}  // namespace ksp_plugin_adapter
}  // namespace principia
