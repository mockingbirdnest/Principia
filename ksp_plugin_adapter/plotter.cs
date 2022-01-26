using System;
using System.Collections.Generic;
using System.Runtime.InteropServices;

namespace principia {
namespace ksp_plugin_adapter {

class Plotter {
  public Plotter(PrincipiaPluginAdapter adapter) {
    adapter_ = adapter;
  }

  public void PlotTrajectories(DisposablePlanetarium planetarium,
                               string main_vessel_guid,
                               double history_length) {
    PlotCelestialTrajectories(planetarium, main_vessel_guid, history_length);
    if (main_vessel_guid == null) {
      return;
    }
    PlotVesselTrajectories(planetarium, main_vessel_guid, history_length);
  }

  private void PlotVesselTrajectories(DisposablePlanetarium planetarium,
                                      string main_vessel_guid,
                                      double history_length) {
    {
      planetarium.PlanetariumPlotPsychohistory(
          Plugin,
          main_vessel_guid,
          history_length,
          VertexBuffer.data,
          VertexBuffer.size,
          out int vertex_count);
      DrawLineMesh(psychohistory_mesh_, vertex_count, adapter_.history_colour,
                   adapter_.history_style);
    }
    {
      planetarium.PlanetariumPlotPrediction(
          Plugin,
          main_vessel_guid,
          VertexBuffer.data,
          VertexBuffer.size,
          out int vertex_count);
      DrawLineMesh(prediction_mesh_, vertex_count, adapter_.prediction_colour,
                   adapter_.prediction_style);
    }

    // Target psychohistory and prediction.
    string target_id = FlightGlobals.fetch.VesselTarget?.GetVessel()?.id.
        ToString();
    if (FlightGlobals.ActiveVessel != null &&
        !adapter_.plotting_frame_selector_.target_frame_selected &&
        target_id != null &&
        Plugin.HasVessel(target_id)) {
      {
        planetarium.PlanetariumPlotPsychohistory(
            Plugin,
            target_id,
            history_length,
            VertexBuffer.data,
            VertexBuffer.size,
            out int vertex_count);
        DrawLineMesh(target_psychohistory_mesh_, vertex_count,
                     adapter_.target_history_colour,
                     adapter_.target_history_style);
      }
      {
        planetarium.PlanetariumPlotPrediction(
            Plugin,
            target_id,
            VertexBuffer.data,
            VertexBuffer.size,
            out int vertex_count);
        DrawLineMesh(target_prediction_mesh_, vertex_count,
                     adapter_.target_prediction_colour,
                     adapter_.target_prediction_style);
      }
    }

    // Main vessel flight plan.
    if (Plugin.FlightPlanExists(main_vessel_guid)) {
      int number_of_segments =
          Plugin.FlightPlanNumberOfSegments(main_vessel_guid);
      for (int i = flight_plan_segment_meshes_.Count; i < number_of_segments; ++i) {
        flight_plan_segment_meshes_.Add(MakeDynamicMesh());
      }
      for (int i = 0; i < number_of_segments; ++i) {
        bool is_burn = i % 2 == 1;
        planetarium.PlanetariumPlotFlightPlanSegment(
            Plugin,
            main_vessel_guid,
            i,
            VertexBuffer.data,
            VertexBuffer.size,
            out int vertex_count);
        DrawLineMesh(flight_plan_segment_meshes_[i],
                     vertex_count,
                     is_burn ? adapter_.burn_colour
                             : adapter_.flight_plan_colour,
                     is_burn ? adapter_.burn_style
                             : adapter_.flight_plan_style);
      }
    }
  }

  private void PlotCelestialTrajectories(DisposablePlanetarium planetarium,
                                         string main_vessel_guid,
                                         double history_length) {
    const double degree = Math.PI / 180;
    UnityEngine.Camera camera = PlanetariumCamera.Camera;
    float vertical_fov = camera.fieldOfView;
    float horizontal_fov =
        UnityEngine.Camera.VerticalToHorizontalFieldOfView(
            vertical_fov, camera.aspect);
    // The angle subtended by the pixel closest to the centre of the viewport.
    double tan_angular_resolution = Math.Min(
            Math.Tan(vertical_fov * degree / 2) / (camera.pixelHeight / 2),
            Math.Tan(horizontal_fov * degree / 2) / (camera.pixelWidth / 2));
    PlotSubtreeTrajectories(planetarium, main_vessel_guid, history_length,
                            Planetarium.fetch.Sun, tan_angular_resolution);
  }

  // Plots the trajectories of |root| and its natural satellites.
  private void PlotSubtreeTrajectories(DisposablePlanetarium planetarium,
                                       string main_vessel_guid,
                                       double history_length,
                                       CelestialBody root,
                                       double tan_angular_resolution) {
    CelestialTrajectories trajectories;
    if (!celestial_trajectory_meshes_.TryGetValue(root, out trajectories)) {
      trajectories = celestial_trajectory_meshes_[root] =
          new CelestialTrajectories();
    }
    var colour = root.orbitDriver?.Renderer?.orbitColor ??
        XKCDColors.SunshineYellow;
    var camera_world_position = ScaledSpace.ScaledToLocalSpace(
        PlanetariumCamera.fetch.transform.position);
    double min_distance_from_camera =
        (root.position - camera_world_position).magnitude;
    if (!adapter_.plotting_frame_selector_.FixesBody(root)) {
      {
        planetarium.PlanetariumPlotCelestialPastTrajectory(
            Plugin,
            root.flightGlobalsIndex,
            history_length,
            VertexBuffer.data,
            VertexBuffer.size,
            out double min_past_distance,
            out int vertex_count);
        min_distance_from_camera =
            Math.Min(min_distance_from_camera, min_past_distance);
        DrawLineMesh(trajectories.past, vertex_count, colour,
                     GLLines.Style.Faded);
      }

      if (main_vessel_guid != null) {
        planetarium.PlanetariumPlotCelestialFutureTrajectory(
            Plugin,
            root.flightGlobalsIndex,
            main_vessel_guid,
            VertexBuffer.data,
            VertexBuffer.size,
            out double min_future_distance,
            out int vertex_count);
        min_distance_from_camera =
            Math.Min(min_distance_from_camera, min_future_distance);
        DrawLineMesh(trajectories.future, vertex_count, colour,
                     GLLines.Style.Solid);
      }
    }
    foreach (CelestialBody child in root.orbitingBodies) {
      // Plot the trajectory of an orbiting body if it could be separated from
      // that of its parent by a pixel of empty space, instead of merely making
      // the line wider.
      if (child.orbit.ApR / min_distance_from_camera >
              2 * tan_angular_resolution) {
        PlotSubtreeTrajectories(planetarium, main_vessel_guid, history_length,
                                child, tan_angular_resolution);
      }
    }
  }

  private void DrawLineMesh(UnityEngine.Mesh mesh,
                            int vertex_count,
                            UnityEngine.Color colour,
                            GLLines.Style style) {
    mesh.vertices = VertexBuffer.vertices;
    var indices = new int[vertex_count];
    for (int i = 0; i < vertex_count; ++i) {
      indices[i] = i;
    }
    var colours = new UnityEngine.Color[VertexBuffer.size];
    if (style == GLLines.Style.Faded) {
      for (int i = 0; i < vertex_count; ++i) {
        var faded_colour = colour;
        // Fade from the opacity of |colour| (when i = 0) down to 20% of that
        // opacity.
        faded_colour.a *= 1 - 0.8f * (i / (float)vertex_count);
        colours[i] = faded_colour;
      }
    } else {
      for (int i = 0; i < vertex_count; ++i) {
        colours[i] = colour;
      }
    }
    mesh.colors = colours;
    mesh.SetIndices(
        indices,
        style == GLLines.Style.Dashed ? UnityEngine.MeshTopology.Lines
                                      : UnityEngine.MeshTopology.LineStrip,
        submesh: 0);
    mesh.RecalculateBounds();
    // If the lines are drawn in layer 31 (Vectors), which sounds more
    // appropriate, they vanish when zoomed out.  Layer 9 works; pay no
    // attention to its name.
    UnityEngine.Graphics.DrawMesh(
        mesh,
        UnityEngine.Vector3.zero,
        UnityEngine.Quaternion.identity,
        GLLines.line_material,
        (int)PrincipiaPluginAdapter.UnityLayers.Atmosphere,
        PlanetariumCamera.Camera);
  }

  private static UnityEngine.Mesh MakeDynamicMesh() {
    var result = new UnityEngine.Mesh();
    result.MarkDynamic();
    return result;
  }

  private readonly PrincipiaPluginAdapter adapter_;

  private IntPtr Plugin => adapter_.Plugin();

  private static class VertexBuffer {
    public static IntPtr data => handle_.AddrOfPinnedObject();
    public static int size => vertices_.Length;

    public static UnityEngine.Vector3[] vertices => vertices_;

    private static UnityEngine.Vector3[] vertices_ =
        new UnityEngine.Vector3[10_000];
    private static GCHandle handle_ =
        GCHandle.Alloc(vertices_, GCHandleType.Pinned);
  }

  private class CelestialTrajectories {
    public UnityEngine.Mesh future = MakeDynamicMesh();
    public UnityEngine.Mesh past = MakeDynamicMesh();
  }

  private Dictionary<CelestialBody, CelestialTrajectories>
      celestial_trajectory_meshes_ =
      new Dictionary<CelestialBody, CelestialTrajectories>();
  private UnityEngine.Mesh psychohistory_mesh_ = MakeDynamicMesh();
  private UnityEngine.Mesh prediction_mesh_ = MakeDynamicMesh();
  private List<UnityEngine.Mesh> flight_plan_segment_meshes_ =
      new List<UnityEngine.Mesh>();
  private UnityEngine.Mesh target_psychohistory_mesh_ = MakeDynamicMesh();
  private UnityEngine.Mesh target_prediction_mesh_ = MakeDynamicMesh();
}

}  // namespace ksp_plugin_adapter
}  // namespace principia
