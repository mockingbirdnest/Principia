using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace principia {
namespace ksp_plugin_adapter {

// Immediately after |PlanetariumCamera| (300).  In particular, before the next
// one, |ScaledSpace| (500), which moves scaled space so that the camera is at
// the origin.  If the adjustment is not performed before |ScaledSpace| runs,
// EnvironmentalVisualEnhancements will draw clouds on an incorrect hemisphere
// (not the camera-facing one).
// Note that KSP provides no timing between 300 and 500, so we have to do this
// ourselves.
[UnityEngine.DefaultExecutionOrder(301)]
public class PlanetariumCameraAdjuster : UnityEngine.MonoBehaviour {
  private void LateUpdate() {
    if (!adapter.PluginRunning() || FlightDriver.Pause) {
      return;
    }
    const double degree = Math.PI / 180;
    var reference_rotation =
        (UnityEngine.QuaternionD)adapter.Plugin().CameraReferenceRotation();
    var camera_roll = UnityEngine.QuaternionD.AngleAxis(camera_roll_ / degree,
                                                        Vector3d.forward);
    if (should_transfer_camera_coordinates) {
      UnityEngine.QuaternionD previous_referenced_pivot =
          previous_camera_reference_rotation_ *
          (UnityEngine.QuaternionD)PlanetariumCamera.fetch.GetPivot().rotation *
          camera_roll;
      // Note that we use a single-precision quaternion here because the
      // double-precision one that comes with KSP does not implement Euler
      // angles.
      UnityEngine.Quaternion new_dereferenced_pivot =
          UnityEngine.QuaternionD.Inverse(reference_rotation) *
          previous_referenced_pivot;
      double new_heading = (new_dereferenced_pivot.eulerAngles.y -
                            Planetarium.InverseRotAngle) * degree;
      double new_pitch = new_dereferenced_pivot.eulerAngles.x * degree;
      // The camera cannot be given nonzero roll by camera controls, but we
      // smoothly bring its roll to 0 over the course of a few frames to make
      // the change in orientation continuous, and thus easier to follow:
      // instantly flipping can be confusing if it is a large change, e.g. Earth
      // equator to Uranus equator.
      camera_roll_ =
          ((new_dereferenced_pivot.eulerAngles.z + 180) % 360 - 180) * degree;
      // Unity has a very mad Euler angle convention that has pitch in
      // [0, π/2] ∪ [3π/2, 2π].
      if (new_pitch > Math.PI) {
        new_pitch -= 2 * Math.PI;
      }
      PlanetariumCamera.fetch.camHdg = (float)new_heading;
      PlanetariumCamera.fetch.camPitch = (float)new_pitch;
      // Use the old reference rotation for this frame: since the change to
      // camera heading and pitch has yet to take effect, the new one would
      // result in a weird orientation for one frame.  Similarly, we keep the
      // existing |camera_roll|.
      reference_rotation = previous_camera_reference_rotation_;
      should_transfer_camera_coordinates = false;
    }
    previous_camera_reference_rotation_ = reference_rotation;
    // Both the scaled space and galaxy cameras are used in the flight scene as
    // well as map view; they should not be reoriented there.
    if (MapView.MapIsEnabled) {
      PlanetariumCamera.fetch.GetPivot().rotation =
          reference_rotation *
          (UnityEngine.QuaternionD)PlanetariumCamera.fetch.GetPivot().rotation *
          camera_roll;
      ScaledCamera.Instance.galaxyCamera.transform.rotation =
          reference_rotation *
          (UnityEngine.QuaternionD)ScaledCamera.Instance.galaxyCamera.transform.rotation *
          camera_roll;
    }
    if (camera_roll_ != 0) {
      // TODO(egg): Should we be doing this in LateUpdate?
      const double roll_change_per_frame = 0.1 /*radians*/;
      if (Math.Abs(camera_roll_) < roll_change_per_frame) {
        camera_roll_ = 0;
      } else {
        camera_roll_ += camera_roll_ > 0 ? -roll_change_per_frame
                                         : roll_change_per_frame;
      }
    }
  }

  // True if a discontinuous change to the camera reference rotation occurred
  // (due to a change in plotting frame), in which case we need special handling
  // to keep the motion of the camera continuous.
  public bool should_transfer_camera_coordinates { private get; set; }

  // The camera reference rotation applied during the previous frame; this is
  // used when transfering camera coordinates to preserve continuity.
  private UnityEngine.QuaternionD previous_camera_reference_rotation_;
  // The roll of the camera; this is set when transferring camera coordinates to
  // preserve continuity of orientation, and gradually brought down to 0 so that
  // the camera is horizontal in the new reference frame.
  private double camera_roll_ = 0;

  public PrincipiaPluginAdapter adapter { private get; set; }
}

}  // namespace ksp_plugin_adapter
}  // namespace principia
