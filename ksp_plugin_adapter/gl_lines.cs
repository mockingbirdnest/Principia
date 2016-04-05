using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace principia {
namespace ksp_plugin_adapter {

internal static class GLLines {
  public enum Style {
    SOLID,
    DASHED,
    FADED,
  }

  public static void Draw(Action line_vertices) {
    try {
      UnityEngine.GL.PushMatrix();
      line_material.SetPass(0);
      UnityEngine.GL.LoadPixelMatrix();
      UnityEngine.GL.Begin(UnityEngine.GL.LINES);
      rendering_lines_ = true;

      Vector3d camera = ScaledSpace.ScaledToLocalSpace(
          PlanetariumCamera.Camera.transform.position);
      // Only consider bodies with an angular radius greater than arcsin 1e-3.
      // From the Earth, this would consider the Moon and the Sun, but ignore
      // Jupiter.  The map view camera is wide-angle, so this is probably
      // overkill.
      // In any case we just want to do that in native code reasonably soon, so
      // this does the trick for now.
      hiding_bodies_ =
        (from body in FlightGlobals.Bodies
         where body.Radius * body.Radius >
               (body.position - camera).sqrMagnitude * 1e-6
         select body).ToArray();

      line_vertices();

      hiding_bodies_ = null;

      rendering_lines_ = false;
      UnityEngine.GL.End();
      UnityEngine.GL.PopMatrix();
    } catch (Exception e) {
      Log.Fatal("Exception while drawing lines: " + e.ToString());
    }
  }

  private static bool IsHidden(Vector3d point) {
    Vector3d camera = ScaledSpace.ScaledToLocalSpace(
        PlanetariumCamera.Camera.transform.position);
    foreach (CelestialBody body in hiding_bodies_) {
      Vector3d camera_to_point = point - camera;
      Vector3d camera_to_body = body.position - camera;
      double inner_product = Vector3d.Dot(camera_to_point, camera_to_body);
      double r_squared = body.Radius * body.Radius;
      // The projections on the camera-body axis of |point| and of the horizon
      // have lengths |inner_product| / d and d - r^2/d, where d is the distance
      // between the camera and the body and r is the body's radius, thus if
      // |inner_product| < d^2 - r^2, |point| is above the plane passing
      // through the horizon.
      // Otherwise, we check whether |point| is within the cone hidden from the
      // camera, by comparing the squared cosines multiplied by
      // d^2|camera_to_point|^2.
      // In addition, we check whether we're inside the body (this covers the
      // cap above the horizon plane and below the surface of the body, which
      // would otherwise be displayed).
      double d_squared_minus_r_squared =
          camera_to_body.sqrMagnitude - r_squared;
      if ((body.position - point).sqrMagnitude < r_squared ||
          (inner_product > d_squared_minus_r_squared &&
           inner_product * inner_product >
               camera_to_point.sqrMagnitude * d_squared_minus_r_squared)) {
        return true;
      }
    }
    return false;
  }

  public static void AddSegment(Vector3d world_begin,
                                Vector3d world_end,
                                bool hide_behind_bodies) {
    if (!rendering_lines_) {
      Log.Fatal("|AddSegment| outside of |DrawLines|");
    }
    if (hide_behind_bodies && (IsHidden(world_begin) || IsHidden(world_end))) {
      return;
    }
    var begin = WorldToMapScreen(world_begin);
    var end = WorldToMapScreen(world_end);
    if (begin.z > 0 && end.z > 0) {
      UnityEngine.GL.Vertex3(begin.x, begin.y, 0);
      UnityEngine.GL.Vertex3(end.x, end.y, 0);
    }
  }

  public static void RenderAndDeleteTrajectory(IntPtr trajectory_iterator,
                                               UnityEngine.Color colour,
                                               Style style) {
    try {
      Vector3d? previous_point = null;

      UnityEngine.GL.Color(colour);
      int size = trajectory_iterator.IteratorSize();

      for (int i = 0;
           !trajectory_iterator.IteratorAtEnd();
           trajectory_iterator.IteratorIncrement(), ++i) {
        Vector3d current_point = (Vector3d)trajectory_iterator.IteratorGetXYZ();
        if (previous_point.HasValue) {
          if (style == Style.FADED) {
            colour.a = (float)(4 * i + size) / (float)(5 * size);
            UnityEngine.GL.Color(colour);
          }
          if (style != Style.DASHED || i % 2 == 0) {
            AddSegment(previous_point.Value,
                       current_point,
                       hide_behind_bodies : true);
          }
        }
        previous_point = current_point;
      }
    } finally {
      Interface.IteratorDelete(ref trajectory_iterator);
    }
  }

  private static UnityEngine.Vector3 WorldToMapScreen(Vector3d world) {
    return PlanetariumCamera.Camera.WorldToScreenPoint(
               ScaledSpace.LocalToScaledSpace(world));
  }

  private static bool rendering_lines_ = false;
  private static CelestialBody[] hiding_bodies_;
  private static UnityEngine.Material line_material_;
  private static UnityEngine.Material line_material {
    get {
      if (line_material_ == null) {
        line_material_ = new UnityEngine.Material(
            UnityEngine.Shader.Find("Particles/Additive"));
      }
      return line_material_;
    }
  }
}

}  // namespace ksp_plugin_adapter
}  // namespace principia
