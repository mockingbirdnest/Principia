using System;
using System.Collections.Generic;
using System.Runtime.InteropServices;

namespace principia {
namespace ksp_plugin_adapter {

class Plotter {
  public Plotter(PrincipiaPluginAdapter adapter) {
    adapter_ = adapter;
  }

  public void Plot(DisposablePlanetarium planetarium,
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
    planetarium.PlanetariumGetPsychohistoryVertices(
        Plugin,
        main_vessel_guid,
        history_length,
        VertexBuffer.data,
        VertexBuffer.size,
        out vertex_count_);
    DrawLineMesh(psychohistory_, adapter_.history_colour, adapter_.history_style);
    planetarium.PlanetariumGetPredictionVertices(
        Plugin,
        main_vessel_guid,
        VertexBuffer.data,
        VertexBuffer.size,
        out vertex_count_);
    DrawLineMesh(prediction_, adapter_.prediction_colour, adapter_.prediction_style);

    // Target psychohistory and prediction.
    string target_id = FlightGlobals.fetch.VesselTarget?.GetVessel()?.id.
        ToString();
    if (FlightGlobals.ActiveVessel != null &&
        !adapter_.plotting_frame_selector_.target_frame_selected &&
        target_id != null &&
        Plugin.HasVessel(target_id)) {
      planetarium.PlanetariumGetPsychohistoryVertices(
          Plugin,
          target_id,
          history_length,
          VertexBuffer.data,
          VertexBuffer.size,
          out vertex_count_);
      DrawLineMesh(target_psychohistory_, adapter_.target_history_colour,
                   adapter_.target_history_style);
      planetarium.PlanetariumGetPredictionVertices(
          Plugin,
          target_id,
          VertexBuffer.data,
          VertexBuffer.size,
          out vertex_count_);
      DrawLineMesh(target_prediction_, adapter_.target_prediction_colour,
                   adapter_.target_prediction_style);
    }
    // Main vessel flight plan.
    if (Plugin.FlightPlanExists(main_vessel_guid)) {
      int number_of_segments =
          Plugin.FlightPlanNumberOfSegments(main_vessel_guid);
      for (int i = flight_plan_segments_.Count; i < number_of_segments; ++i) {
        flight_plan_segments_.Add(MakeDynamicMesh());
      }
      for (int i = 0; i < number_of_segments; ++i) {
        bool is_burn = i % 2 == 1;
        planetarium.PlanetariumGetFlightPlanSegmentVertices(
            Plugin,
            main_vessel_guid,
            i,
            VertexBuffer.data,
            VertexBuffer.size,
            out vertex_count_);
        DrawLineMesh(flight_plan_segments_[i],
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

  private void PlotSubtreeTrajectories(DisposablePlanetarium planetarium,
                                             string main_vessel_guid,
                                             double history_length,
                                             CelestialBody root,
                                             double tan_angular_resolution) {
    CelestialTrajectories trajectories;
    if (!celestial_trajectories_.TryGetValue(root, out trajectories)) {
      trajectories = celestial_trajectories_[root] =
          new CelestialTrajectories();
    }
    var colour = root.orbitDriver?.Renderer?.orbitColor ??
        XKCDColors.SunshineYellow;
    var camera_world_position = ScaledSpace.ScaledToLocalSpace(
        PlanetariumCamera.fetch.transform.position);
    double min_distance_from_camera =
        (root.position - camera_world_position).magnitude;
    if (!adapter_.plotting_frame_selector_.FixesBody(root)) {
      planetarium.PlanetariumGetCelestialPastTrajectoryVertices(
          Plugin,
          root.flightGlobalsIndex,
          history_length,
          VertexBuffer.data,
          VertexBuffer.size,
          out double min_past_distance,
          out vertex_count_);
      min_distance_from_camera =
          Math.Min(min_distance_from_camera, min_past_distance);
      DrawLineMesh(trajectories.past, colour, GLLines.Style.Faded);
      if (main_vessel_guid != null) {
        planetarium.PlanetariumGetCelestialFutureTrajectoryVertices(
            Plugin,
            root.flightGlobalsIndex,
            main_vessel_guid,
            VertexBuffer.data,
            VertexBuffer.size,
            out double min_future_distance,
            out vertex_count_);
        min_distance_from_camera =
            Math.Min(min_distance_from_camera, min_future_distance);
        DrawLineMesh(trajectories.future, colour, GLLines.Style.Solid);
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
                            UnityEngine.Color colour,
                            GLLines.Style style) {
    mesh.vertices = VertexBuffer.get;
    var indices = new int[vertex_count_];
    for (int i = 0; i < vertex_count_; ++i) {
      indices[i] = i;
    }
    var colours = new UnityEngine.Color[VertexBuffer.size];
    if (style == GLLines.Style.Faded) {
      for (int i = 0; i < vertex_count_; ++i) {
        var faded_colour = colour;
        // Fade from the opacity of |colour| (when i = 0) down to 1/4 of that
        // opacity.
        faded_colour.a *= 1 - (float)(4 * i) / (float)(5 * vertex_count_);
        colours[i] = faded_colour;
      }
    } else {
      for (int i = 0; i < vertex_count_; ++i) {
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

    public static UnityEngine.Vector3[] get => vertices_;

    private static UnityEngine.Vector3[] vertices_ =
        new UnityEngine.Vector3[10_000];
    private static GCHandle handle_ =
        GCHandle.Alloc(vertices_, GCHandleType.Pinned);
  }

  private int vertex_count_;

  private class CelestialTrajectories {
    public UnityEngine.Mesh future = MakeDynamicMesh();
    public UnityEngine.Mesh past = MakeDynamicMesh();
  }

  private Dictionary<CelestialBody, CelestialTrajectories> celestial_trajectories_ =
      new Dictionary<CelestialBody, CelestialTrajectories>();
  private UnityEngine.Mesh psychohistory_ = MakeDynamicMesh();
  private UnityEngine.Mesh prediction_ = MakeDynamicMesh();
  private List<UnityEngine.Mesh> flight_plan_segments_ = new List<UnityEngine.Mesh>();
  private UnityEngine.Mesh target_psychohistory_ = MakeDynamicMesh();
  private UnityEngine.Mesh target_prediction_ = MakeDynamicMesh();
}

}  // namespace ksp_plugin_adapter
}  // namespace principia
