using System;
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

internal partial struct Interval {
  public static readonly Interval Empty =
      new Interval
          { min = double.PositiveInfinity, max = double.NegativeInfinity };
  public double measure => min > max ? 0 : max - min;
  public double midpoint => min > max ? double.NaN : min + measure / 2;

  public void Include(double value) {
    if (value < min) {
      min = value;
    }
    if (value > max) {
      max = value;
    }
  }

  public bool Contains(double value) {
    return value >= min && value <= max;
  }

  public Interval IntersectedWith(Interval other) {
    return new Interval
        { min = Math.Max(min, other.min), max = Math.Min(max, other.max) };
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
  ulong CentreCity { get; set; }
  ulong[] PrimaryCities { get; set; }
  ulong[] SecondaryCities { get; set; }
}

internal partial class PlottingFrameParameters : IReferenceFrameParameters,
                                                 IConfigNode {
  public PlottingFrameParameters() {
    primary_city = new ulong[]{};
    secondary_city = new ulong[]{};
  }

  public static explicit operator NavigationFrameParameters(
      PlottingFrameParameters p) {
    if ((FrameType)p.extension == FrameType.ROTATING_PULSATING) {
      return null;
    } else {
      return new NavigationFrameParameters{
          extension = p.extension,
          centre_city= p.centre_city,
          PrimaryCities = p.PrimaryCities,
          SecondaryCities = p.SecondaryCities
      };
    }
  }

  public override bool Equals(object obj) {
    return obj is IReferenceFrameParameters other && this == other;
  }

  public override int GetHashCode() =>
      (extension, centre_city, primary_city.DefaultIfEmpty(no_city).First(),
       secondary_city.DefaultIfEmpty(no_city).First()).GetHashCode();

  public void Load(ConfigNode node) {
    if (!Enum.TryParse(node.GetUniqueValue("extension"),
                       out FrameType extension)) {
      Log.Fatal("Bad FrameType " + node.GetUniqueValue("extension"));
    }
    Extension = extension;
    if (!ulong.TryParse(node.GetUniqueValue("centre_city"), out centre_city)) {
      Log.Fatal("Bad centre_city " + node.GetUniqueValue("centre_city"));
    }
    var primary_cities = new List<ulong>();
    foreach (string value in node.GetValues("primary_cities")) {
      if (!ulong.TryParse(value, out ulong i)) {
        Log.Fatal("Bad primary_cities " + value);
      }
      primary_cities.Add(i);
    }
    primary_city = primary_cities.ToArray();
    var secondary_cities = new List<ulong>();
    foreach (string value in node.GetValues("secondary_cities")) {
      if (!ulong.TryParse(value, out ulong i)) {
        Log.Fatal("Bad secondary_cities " + value);
      }
      secondary_cities.Add(i);
    }
    secondary_city = secondary_cities.ToArray();
  }

  public void Save(ConfigNode node) {
    node.AddValue("extension", Extension.ToString());
    node.AddValue("centre_city", centre_city);
    foreach (ulong i in primary_city) {
      node.AddValue("primary_cities", i);
    }
    foreach (ulong i in secondary_city) {
      node.AddValue("secondary_cities", i);
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
             left.centre_city == right.CentreCity &&
             left.PrimaryCities.SequenceEqual(right.PrimaryCities) &&
             left.SecondaryCities.SequenceEqual(right.SecondaryCities);
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

  public ulong CentreCity {
    get => centre_city;
    set => centre_city = value;
  }

  public ulong[] PrimaryCities {
    get => primary_city;
    set => primary_city = value;
  }

  public ulong[] SecondaryCities {
    get => secondary_city;
    set => secondary_city = value;
  }

  private const ulong no_city = 0xFFFFFFFFFFFFFFFF;
}

internal partial class NavigationFrameParameters : IReferenceFrameParameters {
  public static implicit operator PlottingFrameParameters(
      NavigationFrameParameters p) {
    return new PlottingFrameParameters{
        extension = p.extension,
        centre_city = p.centre_city,
        PrimaryCities = p.PrimaryCities,
        SecondaryCities = p.SecondaryCities
    };
  }

  public override bool Equals(object obj) {
    return obj is IReferenceFrameParameters other && this == other;
  }

  public override int GetHashCode() =>
      (extension, centre_city, primary_city, secondary_city).GetHashCode();

  public static bool operator ==(NavigationFrameParameters left,
                                 IReferenceFrameParameters right) {
    if ((object)left == null && (object)right == null) {
      return true;
    } else if ((object)left == null || (object)right == null) {
      return false;
    } else {
      return left.Extension == right.Extension &&
             left.CentreCity == right.CentreCity &&
             left.PrimaryCities.SequenceEqual(right.PrimaryCities) &&
             left.SecondaryCities.SequenceEqual(right.SecondaryCities);
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

  public ulong CentreCity {
    get => centre_city;
    set => centre_city = value;
  }

  public ulong[] PrimaryCities {
    get => primary_city == no_city ? [] : [ primary_city ];
    set => primary_city = value.DefaultIfEmpty(no_city).Single();
  }

  public ulong[] SecondaryCities {
    get => secondary_city == no_city ? [] : [ secondary_city ];
    set => secondary_city = value.DefaultIfEmpty(no_city).Single();
  }

  private const ulong no_city = 0xFFFFFFFFFFFFFFFF;
}

internal static partial class Interface {
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

  private static Symbols symbols_ = null;

  internal static void LoadSymbols() {
    symbols_ = new Symbols();
  }

  // UnityEngine.Color32 is a union with an int rgba, but it is private, so we
  // use our own union.
  [StructLayout(LayoutKind.Explicit)]
  private struct RGBA32 {
    [FieldOffset(0)]
    public uint rgba;
    [FieldOffset(0)]
    public UnityEngine.Color32 color;
  }

  internal static uint rgba(this UnityEngine.Color color) {
    return new RGBA32{color = color}.rgba;
  }

  internal static UnityEngine.Color32 FromRGBA(uint color) {
    return new RGBA32{rgba = color}.color;
  }
}

} // namespace ksp_plugin_adapter
} // namespace principia
