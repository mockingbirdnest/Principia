using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace principia {
namespace ksp_plugin_adapter {

class BurnEditor {
  public BurnEditor(PrincipiaPluginAdapter adapter,
                    IntPtr plugin,
                    Vessel vessel,
                    double initial_time) {
    Δv_tangent_ = new DifferentialSlider(label            : "Δv tangent",
                                         unit             : "m / s",
                                         log10_lower_rate : Log10ΔvLowerRate,
                                         log10_upper_rate : Log10ΔvUpperRate);
    Δv_normal_ = new DifferentialSlider(label            : "Δv normal",
                                        unit             : "m / s",
                                        log10_lower_rate : Log10ΔvLowerRate,
                                        log10_upper_rate : Log10ΔvUpperRate);
    Δv_binormal_ = new DifferentialSlider(label            : "Δv binormal",
                                          unit             : "m / s",
                                          log10_lower_rate : Log10ΔvLowerRate,
                                          log10_upper_rate : Log10ΔvUpperRate);
    initial_time_ =
        new DifferentialSlider(
                label            : "t initial",
                unit             : null,
                log10_lower_rate : Log10TimeLowerRate,
                log10_upper_rate : Log10TimeUpperRate,
                min_value        : 0,
                max_value        : double.PositiveInfinity,
                formatter        : value =>
                    FlightPlanner.FormatTimeSpan(
                        TimeSpan.FromSeconds(
                            Planetarium.GetUniversalTime() - value)));
    initial_time_.value = initial_time;
    reference_frame_selector_ = new ReferenceFrameSelector(
                                    adapter,
                                    plugin,
                                    ReferenceFrameChanged,
                                    "Manœuvring frame");
    plugin_ = plugin;
    vessel_ = vessel;
    adapter_ = adapter;
    reference_frame_selector_.Reset(
        adapter_.plotting_frame_selector_.get().FrameParameters());
    ComputeEngineCharacteristics();
  }

  // Renders the |BurnEditor|.  Returns true if and only if the settings were
  // changed.
  public bool Render(bool enabled) {
    var old_skin = UnityEngine.GUI.skin;
    UnityEngine.GUI.skin = null;
    UnityEngine.GUILayout.BeginVertical();
    bool changed = false;
    if (enabled) {
      UnityEngine.GUILayout.BeginHorizontal();
      if (UnityEngine.GUILayout.Button("Active Engines")) {
        engine_warning_ = "";
        ComputeEngineCharacteristics();
        changed = true;
      } else if (UnityEngine.GUILayout.Button("Active RCS")) {
        engine_warning_ = "";
        ComputeRCSCharacteristics();
        changed = true;
      } else if (UnityEngine.GUILayout.Button("Instant Impulse")) {
        engine_warning_ = "";
        UseTheForceLuke();
        changed = true;
      }
      UnityEngine.GUILayout.EndHorizontal();
      UnityEngine.GUILayout.TextArea(engine_warning_);
      reference_frame_selector_.RenderButton();
    } else {
      reference_frame_selector_.Hide();
    }
    string frame_warning = "";
    if (!reference_frame_selector_.FrameParameters().Equals(
            adapter_.plotting_frame_selector_.get().FrameParameters())) {
      frame_warning = "Manœuvre frame differs from plotting frame";
    }
    UnityEngine.GUILayout.TextArea(frame_warning);
    changed |= Δv_tangent_.Render(enabled);
    changed |= Δv_normal_.Render(enabled);
    changed |= Δv_binormal_.Render(enabled);
    changed |= initial_time_.Render(enabled);
    changed |= changed_reference_frame_;
    UnityEngine.GUILayout.Label("Burn duration : " + Math.Round(duration_) +
                                " s");
    changed_reference_frame_ = false;
    UnityEngine.GUILayout.EndVertical();
    UnityEngine.GUI.skin = old_skin;
    return changed && enabled;
  }

  public void Reset(NavigationManoeuvre manoeuvre) {
    Burn burn = manoeuvre.burn;
    Δv_tangent_.value = burn.delta_v.x;
    Δv_normal_.value = burn.delta_v.y;
    Δv_binormal_.value = burn.delta_v.z;
    initial_time_.value = burn.initial_time;
    reference_frame_selector_.Reset(burn.frame);
    duration_ = manoeuvre.duration;
  }

  public Burn Burn() {
    return new Burn{
        thrust_in_kilonewtons = thrust_in_kilonewtons_,
        specific_impulse_in_seconds_g0 = specific_impulse_in_seconds_g0_,
        frame = reference_frame_selector_.FrameParameters(),
        initial_time = initial_time_.value,
        delta_v = new XYZ{x = Δv_tangent_.value,
                          y = Δv_normal_.value,
                          z = Δv_binormal_.value}};
  }

  public void ReferenceFrameChanged(NavigationFrameParameters parameters) {
    changed_reference_frame_ = true;
  }

  public void Close() {
    reference_frame_selector_.Dispose();
  }

  private void ComputeEngineCharacteristics() {
    ModuleEngines[] active_engines =
        (from part in vessel_.parts
         select (from PartModule module in part.Modules
                 where module is ModuleEngines &&
                       (module as ModuleEngines).EngineIgnited
                 select module as ModuleEngines)).SelectMany(x => x).ToArray();
    Vector3d reference_direction = vessel_.ReferenceTransform.up;
    double[] thrusts =
        (from engine in active_engines
         select engine.maxThrust *
             (from transform in engine.thrustTransforms
              select Math.Max(0,
                              Vector3d.Dot(reference_direction,
                                           -transform.forward))).Average()).
            ToArray();
    thrust_in_kilonewtons_ = thrusts.Sum();

    // This would use zip if we had 4.0 or later.  We loop for now.
    double Σ_f_over_i_sp = 0;
    for (int i = 0; i < active_engines.Count(); ++i) {
      Σ_f_over_i_sp +=
          thrusts[i] / active_engines[i].atmosphereCurve.Evaluate(0);
    }
    specific_impulse_in_seconds_g0_ = thrust_in_kilonewtons_ / Σ_f_over_i_sp;

    // If there are no engines, fall back onto RCS.
    if (thrust_in_kilonewtons_ == 0) {
      engine_warning_ = "No active engines, falling back to RCS";
      ComputeRCSCharacteristics();
    }
  }

  private void ComputeRCSCharacteristics() {
    ModuleRCS[] active_rcs =
        (from part in vessel_.parts
         select (from PartModule module in part.Modules
                 where module is ModuleRCS &&
                       (module as ModuleRCS).rcsEnabled
                 select module as ModuleRCS)).SelectMany(x => x).ToArray();
    Vector3d reference_direction = vessel_.ReferenceTransform.up;
    // NOTE(egg): NathanKell informs me that in >= 1.0.5, RCS has a useZaxis
    // property, that controls whether they thrust -up or -forward.  The madness
    // keeps piling up.
    double[] thrusts =
        (from engine in active_rcs
         select engine.thrusterPower *
             (from transform in engine.thrusterTransforms
              select Math.Max(0,
                              Vector3d.Dot(reference_direction,
                                           -transform.up))).Average()).
            ToArray();
    thrust_in_kilonewtons_ = thrusts.Sum();

    // This would use zip if we had 4.0 or later.  We loop for now.
    double Σ_f_over_i_sp = 0;
    for (int i = 0; i < active_rcs.Count(); ++i) {
      Σ_f_over_i_sp +=
          thrusts[i] / active_rcs[i].atmosphereCurve.Evaluate(0);
    }
    specific_impulse_in_seconds_g0_ = thrust_in_kilonewtons_ / Σ_f_over_i_sp;

    // If RCS provides no thrust, model a virtually instant burn.
    if (thrust_in_kilonewtons_ == 0) {
      engine_warning_ = "No active RCS, modeling as instant burn";
      UseTheForceLuke();
    }
  }

  private void UseTheForceLuke() {
    thrust_in_kilonewtons_ = 1E15;
    specific_impulse_in_seconds_g0_ = 3E7;
  }

  private DifferentialSlider Δv_tangent_;
  private DifferentialSlider Δv_normal_;
  private DifferentialSlider Δv_binormal_;
  private DifferentialSlider initial_time_;
  private ReferenceFrameSelector reference_frame_selector_;
  private double thrust_in_kilonewtons_;
  private double specific_impulse_in_seconds_g0_;
  private double duration_;

  private const double Log10ΔvLowerRate = -3.0;
  private const double Log10ΔvUpperRate = 3.5;
  private const double Log10TimeLowerRate = 0.0;
  private const double Log10TimeUpperRate = 7.0;

  // Not owned.
  private readonly IntPtr plugin_;
  private readonly Vessel vessel_;
  private readonly PrincipiaPluginAdapter adapter_;

  private bool changed_reference_frame_ = false;
  private string engine_warning_ = "";
}

}  // namespace ksp_plugin_adapter
}  // namespace principia
