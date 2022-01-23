using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace principia {
namespace ksp_plugin_adapter {

class Plotter {
  public Plotter(PrincipiaPluginAdapter adapter) {
    adapter_ = adapter;
  }

  public unsafe void Plot(DisposablePlanetarium planetarium,
                          string main_vessel_guid,
                          double history_length) {
    fixed (UnityEngine.Vector3* vertices_data = vertices_) {
      vertices_data_ = (IntPtr)vertices_data;
      PlotCelestialTrajectories(planetarium, main_vessel_guid, history_length);
      if (main_vessel_guid == null) {
        return;
      }
      PlotVesselTrajectories(planetarium, main_vessel_guid, history_length);
      vertices_data_ = IntPtr.Zero;
    }
  }

  private void PlotVesselTrajectories(DisposablePlanetarium planetarium,
                                      string main_vessel_guid,
                                      double history_length) {
    planetarium.PlanetariumGetPsychohistoryVertices(
        Plugin,
        main_vessel_guid,
        history_length,
        vertices_data_,
        vertices_.Length,
        out vertex_count_);
    DrawLineMesh(psychohistory_, adapter_.history_colour, adapter_.history_style);
    planetarium.PlanetariumGetPredictionVertices(
        Plugin,
        main_vessel_guid,
        vertices_data_,
        vertices_.Length,
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
          vertices_data_,
          vertices_.Length,
          out vertex_count_);
      DrawLineMesh(target_psychohistory_, adapter_.target_history_colour,
                   adapter_.target_history_style);
    }
    // Main vessel flight plan.
    if (Plugin.FlightPlanExists(main_vessel_guid)) {
      int number_of_anomalous_manœuvres =
          Plugin.FlightPlanNumberOfAnomalousManoeuvres(main_vessel_guid);
      int number_of_manœuvres =
          Plugin.FlightPlanNumberOfManoeuvres(main_vessel_guid);
      int number_of_segments =
          Plugin.FlightPlanNumberOfSegments(main_vessel_guid);
      for (int i = flight_plan_segments_.Count; i < number_of_segments; ++i) {
        flight_plan_segments_.Add(new UnityEngine.Mesh());
      }
      for (int i = 0; i < number_of_segments; ++i) {
        bool is_burn = i % 2 == 1;
        planetarium.PlanetariumGetFlightPlanSegmentVertices(
            Plugin,
            main_vessel_guid,
            i,
            vertices_data_,
            vertices_.Length,
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
    MainWindow.trace = $"{UnityEngine.Camera.current.name}, {PlanetariumCamera.Camera.name} {ScaledSpace.LocalToScaledSpace(Vector3d.zero)}";
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
          vertices_data_,
          vertices_.Length,
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
            vertices_data_,
            vertices_.Length,
            out double min_future_distance,
            out vertex_count_);
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
    filled_vertices_ |= vertex_count_ == vertices_.Length;
    if (!MainWindow.reuse_meshes) {
      mesh = new UnityEngine.Mesh();
    }
    mesh.vertices = vertices_;
    var indices = new int[vertex_count_];
    for (int i = 0; i < vertex_count_; ++i) {
      indices[i] = i;
    }
    var colours = new UnityEngine.Color[vertices_.Length];
    if (style == GLLines.Style.Faded) {
      for (int i = 0; i < colours.Length; ++i) {
        var faded_colour = colour;
        // Fade from the opacity of |colour| (when i = 0) down to 1/4 of that
        // opacity.
        faded_colour.a *= 1 - (float)(4 * i) / (float)(5 * colours.Length);
        colours[i] = faded_colour;
      }
    } else {
      for (int i = 0; i < colours.Length; ++i) {
        colours[i] = colour;
      }
    }
    mesh.colors = colours;
    mesh.SetIndices(
        indices,
        style == GLLines.Style.Dashed ? UnityEngine.MeshTopology.Lines
                                      : UnityEngine.MeshTopology.LineStrip,
        submesh: 0);
    if (MainWindow.reuse_meshes) {
      mesh.RecalculateBounds();
    }
    if (MainWindow.now) {
      UnityEngine.Graphics.DrawMeshNow(mesh, UnityEngine.Vector3.zero, UnityEngine.Quaternion.identity);
    } else {
      UnityEngine.Graphics.DrawMesh(
          mesh,
          UnityEngine.Vector3.zero,
          UnityEngine.Quaternion.identity,
          GLLines.line_material,
          MainWindow.layer,
          PlanetariumCamera.Camera);
    }
  }

  private readonly PrincipiaPluginAdapter adapter_;

  private IntPtr Plugin => adapter_.Plugin();

  private UnityEngine.Vector3[] vertices_ = new UnityEngine.Vector3[10_000];
  // A pointer to vertices_; vertices_ should be pinned and should not be
  // reassigned while this is non-null.
  private IntPtr vertices_data_ = IntPtr.Zero;
  private int vertex_count_;
  private bool filled_vertices_ = false;

  private class CelestialTrajectories {
    public UnityEngine.Mesh future = new UnityEngine.Mesh();
    public UnityEngine.Mesh past = new UnityEngine.Mesh();
  }

  private Dictionary<CelestialBody, CelestialTrajectories> celestial_trajectories_ =
      new Dictionary<CelestialBody, CelestialTrajectories>();
  private UnityEngine.Mesh psychohistory_ = new UnityEngine.Mesh();
  private UnityEngine.Mesh prediction_ = new UnityEngine.Mesh();
  private List<UnityEngine.Mesh> flight_plan_segments_ = new List<UnityEngine.Mesh>();
  private UnityEngine.Mesh target_psychohistory_ = new UnityEngine.Mesh();
  private UnityEngine.Mesh target_prediction_ = new UnityEngine.Mesh();
}

}  // namespace ksp_plugin_adapter
}  // namespace principia
