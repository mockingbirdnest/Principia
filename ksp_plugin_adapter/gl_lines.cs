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

      Vector3d camera = ScaledSpace.ScaledToLocalSpace(
          PlanetariumCamera.Camera.transform.position);
      line_vertices();

      UnityEngine.GL.End();
      UnityEngine.GL.PopMatrix();
    } catch (Exception e) {
      Log.Fatal("Exception while drawing lines: " + e.ToString());
    }
  }

  public static void AddSegment(Vector3d world_begin,
                                Vector3d world_end) {
    var begin = WorldToMapScreen(world_begin);
    var end = WorldToMapScreen(world_end);
    if (begin.z > 0 && end.z > 0) {
      UnityEngine.GL.Vertex3(begin.x, begin.y, 0);
      UnityEngine.GL.Vertex3(end.x, end.y, 0);
    }
  }

  public static IntPtr NewPlanetarium(IntPtr plugin,
                                       XYZ sun_world_position) {
    UnityEngine.Camera camera = PlanetariumCamera.Camera;
    UnityEngine.Vector3 opengl_camera_x_in_world =
        camera.cameraToWorldMatrix.MultiplyVector(
            new UnityEngine.Vector3(1, 0, 0));
    UnityEngine.Vector3 opengl_camera_y_in_world =
        camera.cameraToWorldMatrix.MultiplyVector(
            new UnityEngine.Vector3(0, 1, 0));
    UnityEngine.Vector3 opengl_camera_z_in_world =
        camera.cameraToWorldMatrix.MultiplyVector(
            new UnityEngine.Vector3(0, 0, 1));
    UnityEngine.Vector3 camera_position_in_world =
        ScaledSpace.ScaledToLocalSpace(camera.transform.position);

    // For explanations regarding the OpenGL projection matrix, see
    // http://www.songho.ca/opengl/gl_projectionmatrix.html.  The on-centre
    // projection matrix has the form:
    //   n / w                0                0                0
    //     0                n / h              0                0
    //     0                  0        (n + f) / (n - f)  2 f n / (n - f)
    //     0                  0               -1                0
    // where n and f are the near- and far-clipping distances, and w and h
    // are the half-width and half-height of the screen seen in the focal plane.
    // n is also the focal distance, but we prefer to make that distance 1 metre
    // to avoid having to rescale the result.  The only actual effect of n is
    // the clipping distance, and in space, no one can hear you clip.
    double m00 = camera.projectionMatrix[0, 0];
    double m11 = camera.projectionMatrix[1, 1];
    double field_of_view = Math.Atan2(Math.Sqrt(m00 * m00 + m11 * m11),
                                      m00 * m11);
    return plugin.PlanetariumCreate(
               sun_world_position,
               (XYZ)(Vector3d)opengl_camera_x_in_world,
               (XYZ)(Vector3d)opengl_camera_y_in_world,
               (XYZ)(Vector3d)opengl_camera_z_in_world,
               (XYZ)(Vector3d)camera_position_in_world,
               /*focal=*/1,
               field_of_view);
  }

  public static void PlotRP2Lines(DisposableIterator rp2_lines_iterator,
                                  UnityEngine.Color colour,
                                  Style style) {
    UnityEngine.GL.Color(colour);

    // First evaluate the total size of the lines.
    int size = 0;
    for (;
         !rp2_lines_iterator.IteratorAtEnd();
         rp2_lines_iterator.IteratorIncrement()) {
      using (DisposableIterator rp2_line_iterator =
                rp2_lines_iterator.IteratorGetRP2LinesIterator()) {
        size += rp2_line_iterator.IteratorSize();
      }
    }

    // Reset the iterator and do the actual plotting.
    rp2_lines_iterator.IteratorReset();
    int index = 0;
    for (;
         !rp2_lines_iterator.IteratorAtEnd();
         rp2_lines_iterator.IteratorIncrement()) {
      using (DisposableIterator rp2_line_iterator =
                rp2_lines_iterator.IteratorGetRP2LinesIterator()) {
        XY? previous_rp2_point = null;
        for (;
             !rp2_line_iterator.IteratorAtEnd();
             rp2_line_iterator.IteratorIncrement()) {
          XY current_rp2_point = ToScreen(
              rp2_line_iterator.IteratorGetRP2LineXY());
          if (previous_rp2_point.HasValue) {
            if (style == Style.FADED) {
              colour.a = 1 - (float)(4 * index) / (float)(5 * size);
              UnityEngine.GL.Color(colour);
            }
            if (style != Style.DASHED || index % 2 == 1) {
              UnityEngine.GL.Vertex3((float)previous_rp2_point.Value.x,
                                      (float)previous_rp2_point.Value.y,
                                      0);
              UnityEngine.GL.Vertex3((float)current_rp2_point.x,
                                      (float)current_rp2_point.y,
                                      0);
            }
          }
          previous_rp2_point = current_rp2_point;
          ++index;
        }
      }
    }
  }

  private static UnityEngine.Vector3 WorldToMapScreen(Vector3d world) {
    return PlanetariumCamera.Camera.WorldToScreenPoint(
               ScaledSpace.LocalToScaledSpace(world));
  }

  private static XY ToScreen(XY rp2_point) {
    UnityEngine.Camera camera = PlanetariumCamera.Camera;
    return new XY{x = (rp2_point.x * camera.projectionMatrix[0, 0] + 1.0) *
                      0.5 * camera.pixelWidth,
                  y = (rp2_point.y * camera.projectionMatrix[1, 1] + 1.0) *
                      0.5 * camera.pixelHeight};
   }

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
