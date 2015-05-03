using System;
using System.Collections.Generic;
using System.Runtime.InteropServices;

namespace principia {
namespace ksp_plugin_adapter {

public partial class PrincipiaPluginAdapter : ScenarioModule {

  internal const string kDllPath = "GameData/Principia/principia.dll";

  [StructLayout(LayoutKind.Sequential)]
  private struct XYZ {
    public double x, y, z;
    public static explicit operator XYZ(Vector3d v) {
      return new XYZ{x = v.x, y = v.y, z = v.z};
    }
    public static explicit operator Vector3d(XYZ v) {
      return new Vector3d{x = v.x, y = v.y, z = v.z};
    }
  };

  [StructLayout(LayoutKind.Sequential)]
  private struct WXYZ {
    public double w, x, y, z;
    public static explicit operator WXYZ(UnityEngine.QuaternionD q) {
      return new WXYZ{w = q.w, x = q.x, y = q.y, z = q.z};
    }
    public static explicit operator UnityEngine.QuaternionD(WXYZ q) {
      return new UnityEngine.QuaternionD{w = q.w, x = q.x, y = q.y, z = q.z};
    }
  };

  [StructLayout(LayoutKind.Sequential)]
  private struct LineSegment {
    public XYZ begin, end;
  };

  [StructLayout(LayoutKind.Sequential)]
  private struct QP {
    public XYZ q, p;
  };

  [StructLayout(LayoutKind.Sequential)]
  struct KSPPart {
    public XYZ world_position;
    public XYZ world_velocity;
    public double mass;
    public XYZ gravitational_acceleration_to_be_applied_by_ksp;
    public uint id;
  };

  // Plugin interface.

  [DllImport(dllName           : kDllPath,
             EntryPoint        = "principia__SayHello",
             CallingConvention = CallingConvention.Cdecl)]
  private static extern IntPtr SayHello();

  [DllImport(dllName           : kDllPath,
             EntryPoint        = "principia__NewPlugin",
             CallingConvention = CallingConvention.Cdecl)]
  private static extern IntPtr NewPlugin(
      double initial_time,
      int sun_index,
      double sun_gravitational_parameter,
      double planetarium_rotation_in_degrees);

  [DllImport(dllName           : kDllPath,
             EntryPoint        = "principia__DeletePlugin",
             CallingConvention = CallingConvention.Cdecl)]
  private static extern void DeletePlugin(ref IntPtr plugin);

  [DllImport(dllName           : kDllPath,
             EntryPoint        = "principia__InsertCelestial",
             CallingConvention = CallingConvention.Cdecl)]
  private static extern void InsertCelestial(
      IntPtr plugin,
      int celestial_index,
      double gravitational_parameter,
      int parent_index,
      QP from_parent);

  [DllImport(dllName           : kDllPath,
             EntryPoint        = "principia__EndInitialization",
             CallingConvention = CallingConvention.Cdecl)]
  private static extern void EndInitialization(IntPtr plugin);

  [DllImport(dllName           : kDllPath,
             EntryPoint        = "principia__UpdateCelestialHierarchy",
             CallingConvention = CallingConvention.Cdecl)]
  private static extern void UpdateCelestialHierarchy(
      IntPtr plugin,
      int celestial_index,
      int parent_index);

  [DllImport(dllName           : kDllPath,
             EntryPoint        = "principia__InsertOrKeepVessel",
             CallingConvention = CallingConvention.Cdecl)]
  private static extern bool InsertOrKeepVessel(
      IntPtr plugin,
      [MarshalAs(UnmanagedType.LPStr)] String vessel_guid,
      int parent_index);

  [DllImport(dllName           : kDllPath,
             EntryPoint        = "principia__SetVesselStateOffset",
             CallingConvention = CallingConvention.Cdecl)]
  private static extern void SetVesselStateOffset(
      IntPtr plugin,
      [MarshalAs(UnmanagedType.LPStr)] String vessel_guid,
      QP from_parent);

  [DllImport(dllName           : kDllPath,
             EntryPoint        = "principia__AdvanceTime",
             CallingConvention = CallingConvention.Cdecl)]
  private static extern void AdvanceTime(
      IntPtr plugin, 
      double t,
      double planetarium_rotation);

  [DllImport(dllName           : kDllPath,
             EntryPoint        = "principia__VesselFromParent",
             CallingConvention = CallingConvention.Cdecl)]
  private static extern QP VesselFromParent(
      IntPtr plugin,
      [MarshalAs(UnmanagedType.LPStr)] String vessel_guid);

  [DllImport(dllName           : kDllPath,
             EntryPoint        = "principia__CelestialFromParent",
             CallingConvention = CallingConvention.Cdecl)]
  private static extern QP CelestialFromParent(
      IntPtr plugin,
      int celestial_index);

  [DllImport(dllName           : kDllPath,
             EntryPoint        =
                 "principia__NewBodyCentredNonRotatingTransforms",
             CallingConvention = CallingConvention.Cdecl)]
  private static extern IntPtr NewBodyCentredNonRotatingTransforms(
      IntPtr plugin,
      int reference_body_index);

  [DllImport(dllName           : kDllPath,
             EntryPoint        = "principia__NewBarycentricRotatingTransforms",
             CallingConvention = CallingConvention.Cdecl)]
  private static extern IntPtr NewBarycentricRotatingTransforms(
      IntPtr plugin,
      int primary_index,
      int secondary_index);

  [DllImport(dllName           : kDllPath,
             EntryPoint        = "principia__DeleteTransforms",
             CallingConvention = CallingConvention.Cdecl)]
  private static extern void DeleteTransforms(ref IntPtr transforms);

  [DllImport(dllName           : kDllPath,
             EntryPoint        = "principia__RenderedVesselTrajectory",
             CallingConvention = CallingConvention.Cdecl)]
  private static extern IntPtr RenderedVesselTrajectory(
      IntPtr plugin,
      [MarshalAs(UnmanagedType.LPStr)] String vessel_guid,
      IntPtr transforms,
      XYZ sun_world_position);

  [DllImport(dllName           : kDllPath,
             EntryPoint        = "principia__RenderedPrediction",
             CallingConvention = CallingConvention.Cdecl)]
  private static extern IntPtr RenderedPrediction(
      IntPtr plugin,
      IntPtr transforms,
      XYZ sun_world_position);

  [DllImport(dllName           : kDllPath,
             EntryPoint        = "principia__NumberOfSegments",
             CallingConvention = CallingConvention.Cdecl)]
  private static extern int NumberOfSegments(IntPtr line);

  [DllImport(dllName           : kDllPath,
             EntryPoint        = "principia__FetchAndIncrement",
             CallingConvention = CallingConvention.Cdecl)]
  private static extern LineSegment FetchAndIncrement(IntPtr line);

  [DllImport(dllName           : kDllPath,
             EntryPoint        = "principia__AtEnd",
             CallingConvention = CallingConvention.Cdecl)]
  private static extern bool AtEnd(IntPtr line);

  [DllImport(dllName           : kDllPath,
             EntryPoint        = "principia__DeleteLineAndIterator",
             CallingConvention = CallingConvention.Cdecl)]
  private static extern void DeleteLineAndIterator(ref IntPtr line);

  [DllImport(dllName           : kDllPath,
             EntryPoint        = "principia__set_predicted_vessel",
             CallingConvention = CallingConvention.Cdecl)]
  private static extern void set_predicted_vessel(
      IntPtr plugin,
      [MarshalAs(UnmanagedType.LPStr)] String vessel_guid);

  [DllImport(dllName           : kDllPath,
             EntryPoint        = "principia__clear_predicted_vessel",
             CallingConvention = CallingConvention.Cdecl)]
  private static extern void clear_predicted_vessel(IntPtr plugin);

  [DllImport(dllName           : kDllPath,
             EntryPoint        = "principia__set_prediction_length",
             CallingConvention = CallingConvention.Cdecl)]
  private static extern void set_prediction_length(IntPtr plugin, double t);

  [DllImport(dllName           : kDllPath,
             EntryPoint        = "principia__set_prediction_step",
             CallingConvention = CallingConvention.Cdecl)]
  private static extern void set_prediction_step(IntPtr plugin, double t);

  [DllImport(dllName             : kDllPath,
             EntryPoint =        "principia__has_vessel",
             CallingConvention = CallingConvention.Cdecl)]
  private static extern bool has_vessel(
      IntPtr plugin,
      [MarshalAs(UnmanagedType.LPStr)] String vessel_guid);

  [DllImport(dllName: kDllPath,
             EntryPoint        = "principia__AddVesselToNextPhysicsBubble",
             CallingConvention = CallingConvention.Cdecl)]
  private static extern void AddVesselToNextPhysicsBubble(
      IntPtr plugin,
      [MarshalAs(UnmanagedType.LPStr)] String vessel_guid,
      KSPPart[] parts,
      int count);

  [DllImport(dllName           : kDllPath,
             EntryPoint        = "principia__PhysicsBubbleIsEmpty",
             CallingConvention = CallingConvention.Cdecl)]
  private static extern bool PhysicsBubbleIsEmpty(IntPtr plugin);

  [DllImport(dllName           : kDllPath,
             EntryPoint        = "principia__BubbleDisplacementCorrection",
             CallingConvention = CallingConvention.Cdecl)]
  private static extern XYZ BubbleDisplacementCorrection(IntPtr plugin,
                                                         XYZ sun_position);

  [DllImport(dllName           : kDllPath,
             EntryPoint        = "principia__BubbleVelocityCorrection",
             CallingConvention = CallingConvention.Cdecl)]
  private static extern XYZ BubbleVelocityCorrection(IntPtr plugin,
                                                     int reference_body_index);

  [DllImport(dllName           : kDllPath,
             EntryPoint        = "principia__NavballOrientation",
             CallingConvention = CallingConvention.Cdecl)]
  private static extern WXYZ NavballOrientation(
      IntPtr plugin,
      IntPtr transforms,
      XYZ sun_world_position,
      XYZ ship_world_position);

  [DllImport(dllName           : kDllPath,
             EntryPoint        = "principia__VesselTangent",
             CallingConvention = CallingConvention.Cdecl)]
  private static extern XYZ VesselTangent(
      IntPtr plugin,
      [MarshalAs(UnmanagedType.LPStr)] String vessel_guid,
      IntPtr transforms);

  [DllImport(dllName           : kDllPath,
             EntryPoint        = "principia__current_time",
             CallingConvention = CallingConvention.Cdecl)]
  private static extern double current_time(IntPtr plugin);

  [DllImport(dllName           : kDllPath,
             EntryPoint        = "principia__SerializePlugin",
             CallingConvention = CallingConvention.Cdecl)]
  private static extern IntPtr SerializePlugin(IntPtr plugin,
                                               ref IntPtr serializer);

  [DllImport(dllName           : kDllPath,
             EntryPoint        = "principia__DeletePluginSerialization",
             CallingConvention = CallingConvention.Cdecl)]
  private static extern void DeletePluginSerialization(
      ref IntPtr serialization);

  [DllImport(dllName           : kDllPath,
             EntryPoint        = "principia__DeserializePlugin",
             CallingConvention = CallingConvention.Cdecl)]
  private static extern void DeserializePlugin(
      [MarshalAs(UnmanagedType.LPStr)] String serialization,
      int serialization_size,
      ref IntPtr deserializer,
      ref IntPtr plugin);

}

}  // namespace ksp_plugin_adapter
}  // namespace principia
