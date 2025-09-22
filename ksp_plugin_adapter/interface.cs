﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.InteropServices;

namespace principia {
namespace ksp_plugin_adapter {

internal partial class Status {
  public static Status OK = new Status{ error = 0 };

  public bool is_aborted() {
    return error == 10;
  }

  public bool is_deadline_exceeded() {
    return error == 4;
  }

  public bool is_failed_precondition() {
    return error == 9;
  }

  public bool is_invalid_argument() {
    return error == 3;
  }

  public bool is_out_of_range() {
    return error == 11;
  }

  public bool is_resource_exhausted() {
    return error == 8;
  }

  public bool is_unavailable() {
    return error == 14;
  }

  public bool ok() {
    return error == 0;
  }

  public void Update(Status s) {
    if (ok() && !s.ok()) {
      this.error = s.error;
      this.message = s.message;
    }
  }
}

public partial struct XYZ {
  public static explicit operator XYZ(UnityEngine.Vector3 v) {
    return new XYZ{ x = v.x, y = v.y, z = v.z };
  }

  public static explicit operator XYZ(Vector3d v) {
    return new XYZ { x = v.x, y = v.y, z = v.z };
  }

  public static explicit operator UnityEngine.Vector3(XYZ v) {
    return new UnityEngine.Vector3((float)v.x, (float)v.y, (float)v.z);
  }

  public static explicit operator Vector3d(XYZ v) {
    return new Vector3d{ x = v.x, y = v.y, z = v.z };
  }
}

internal partial struct WXYZ {
  public static explicit operator WXYZ(UnityEngine.Quaternion q) {
    return new WXYZ{ w = q.w, x = q.x, y = q.y, z = q.z };
  }

  public static explicit operator WXYZ(UnityEngine.QuaternionD q) {
    return new WXYZ{ w = q.w, x = q.x, y = q.y, z = q.z };
  }

  public static explicit operator UnityEngine.QuaternionD(WXYZ q) {
    return new UnityEngine.QuaternionD{ w = q.w, x = q.x, y = q.y, z = q.z };
  }
}

public enum FrameType {
  BARYCENTRIC_ROTATING = 6001,
  BODY_CENTRED_NON_ROTATING = 6000,
  BODY_CENTRED_PARENT_DIRECTION = 6002,
  BODY_SURFACE = 6003,
  ROTATING_PULSATING = 6004,
}

// We imbue PlottingFrameParameters and NavigationFrameParameters with a common
// interface and make the latter behave like a subtype of the former.

internal interface IReferenceFrameParameters {
  FrameType Extension  { get; set; }
  int CentreIndex { get; set; }
  int[] PrimaryIndices { get; set; }
  int[] SecondaryIndices { get; set; }
}

internal partial class PlottingFrameParameters : IReferenceFrameParameters,
                                                 IConfigNode {
  public PlottingFrameParameters() {
    primary_index = new int[]{};
    secondary_index = new int[]{};
  }

  public static explicit operator NavigationFrameParameters(
      PlottingFrameParameters p) {
    if ((FrameType)p.extension == FrameType.ROTATING_PULSATING) {
      return null;
    } else {
      return new NavigationFrameParameters{
          extension = p.extension,
          centre_index = p.centre_index,
          PrimaryIndices = p.PrimaryIndices,
          SecondaryIndices = p.SecondaryIndices
      };
    }
  }

  public override bool Equals(object obj) {
    return obj is IReferenceFrameParameters other && this == other;
  }

  public override int GetHashCode() =>
      (extension, centre_index, primary_index.DefaultIfEmpty(-1).First(),
       secondary_index.DefaultIfEmpty(-1).First()).GetHashCode();

  public void Load(ConfigNode node) {
    if (!Enum.TryParse(node.GetUniqueValue("extension"),
                       out FrameType extension)) {
      Log.Fatal("Bad FrameType " + node.GetUniqueValue("extension"));
    }
    Extension = extension;
    if (!int.TryParse(node.GetUniqueValue("centre_index"), out centre_index)) {
      Log.Fatal("Bad centre_index " + node.GetUniqueValue("centre_index"));
    }
    var primary_indices = new List<int>();
    foreach (string value in node.GetValues("primary_indices")) {
      if (!int.TryParse(value, out int i)) {
        Log.Fatal("Bad primary_indices " + value);
      }
      primary_indices.Add(i);
    }
    primary_index = primary_indices.ToArray();
    var secondary_indices = new List<int>();
    foreach (string value in node.GetValues("secondary_indices")) {
      if (!int.TryParse(value, out int i)) {
        Log.Fatal("Bad secondary_indices " + value);
      }
      secondary_indices.Add(i);
    }
    secondary_index = secondary_indices.ToArray();
  }

  public void Save(ConfigNode node) {
    node.AddValue("extension", Extension.ToString());
    node.AddValue("centre_index", centre_index);
    foreach (int i in primary_index) {
      node.AddValue("primary_indices", i);
    }
    foreach (int i in secondary_index) {
      node.AddValue("secondary_indices", i);
    }
  }

  public static bool operator ==(PlottingFrameParameters left,
                                 IReferenceFrameParameters right) {
    if ((object)left == null && (object)right == null) {
      return true;
    } else if ((object)left == null || (object)right == null) {
      return false;
    } else {
      return left.Extension == right.Extension &&
             left.centre_index == right.CentreIndex &&
             left.PrimaryIndices.SequenceEqual(right.PrimaryIndices) &&
             left.SecondaryIndices.SequenceEqual(right.SecondaryIndices);
    }
  }

  public static bool operator !=(PlottingFrameParameters left,
                                 IReferenceFrameParameters right) {
    return !(left == right);
  }

  public FrameType Extension {
    get => (FrameType)extension;
    set => extension = (int)value;
  }

  public int CentreIndex {
    get => centre_index;
    set => centre_index = value;
  }

  public int[] PrimaryIndices {
    get => primary_index;
    set => primary_index = value;
  }

  public int[] SecondaryIndices {
    get => secondary_index;
    set => secondary_index = value;
  }
}

internal partial class NavigationFrameParameters : IReferenceFrameParameters {
  public static implicit operator PlottingFrameParameters(
      NavigationFrameParameters p) {
    return new PlottingFrameParameters{
        extension = p.extension,
        centre_index = p.centre_index,
        PrimaryIndices = p.PrimaryIndices,
        SecondaryIndices = p.SecondaryIndices
    };
  }

  public override bool Equals(object obj) {
    return obj is IReferenceFrameParameters other && this == other;
  }

  public override int GetHashCode() =>
      (extension, centre_index, primary_index, secondary_index).GetHashCode();

  public static bool operator ==(NavigationFrameParameters left,
                                 IReferenceFrameParameters right) {
    if ((object)left == null && (object)right == null) {
      return true;
    } else if ((object)left == null || (object)right == null) {
      return false;
    } else {
      return left.Extension == right.Extension &&
             left.CentreIndex == right.CentreIndex &&
             left.PrimaryIndices.SequenceEqual(right.PrimaryIndices) &&
             left.SecondaryIndices.SequenceEqual(right.SecondaryIndices);
    }
  }

  public static bool operator !=(NavigationFrameParameters left,
                                 IReferenceFrameParameters right) {
    return !(left == right);
  }

  public FrameType Extension {
    get => (FrameType)extension;
    set => extension = (int)value;
  }

  public int CentreIndex {
    get => centre_index;
    set => centre_index = value;
  }

  public int[] PrimaryIndices {
    get => primary_index == -1 ? new int[] {} : new[] { primary_index };
    set => primary_index = value.DefaultIfEmpty(-1).Single();
  }

  public int[] SecondaryIndices {
    get => secondary_index == -1 ? new int[] {} : new[] { secondary_index };
    set => secondary_index = value.DefaultIfEmpty(-1).Single();
  }
}

internal static partial class Interface {
  internal const string dll_path = "principia";

  internal static KeplerianElements Elements(this Orbit orbit) {
    double mean_motion = 2 * Math.PI / orbit.period;
    return new KeplerianElements{
        eccentricity                           = orbit.eccentricity,
        semimajor_axis                         = double.NaN,
        mean_motion                            = mean_motion,
        inclination_in_degrees                 = orbit.inclination,
        longitude_of_ascending_node_in_degrees = orbit.LAN,
        argument_of_periapsis_in_degrees       = orbit.argumentOfPeriapsis,
        mean_anomaly                           =
            orbit.meanAnomalyAtEpoch - orbit.epoch * mean_motion
    };
  }

  // These two functions are not journalled, so the interface is handwritten on
  // both sides.

  [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
  internal delegate void ActivateRecorderDelegate(bool activate);

  [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
  private delegate void InitGoogleLoggingDelegate();

  private partial class Symbols {
    public ActivateRecorderDelegate ActivateRecorder =
        Loader.LoadFunction<ActivateRecorderDelegate>(
            "principia__ActivateRecorder");
    public InitGoogleLoggingDelegate InitGoogleLogging =
        Loader.LoadFunction<InitGoogleLoggingDelegate>(
            "principia__InitGoogleLogging");
  }

  internal static void ActivateRecorder(bool activate) {
    symbols_.ActivateRecorder(activate);
  }

  internal static void InitGoogleLogging() {
    symbols_.InitGoogleLogging();
  }

  private static Symbols symbols_ = null;

  internal static void LoadSymbols() {
    symbols_ = new Symbols();
  }
}

} // namespace ksp_plugin_adapter
} // namespace principia
