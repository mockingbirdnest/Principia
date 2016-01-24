using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace principia {
namespace ksp_plugin_adapter {

class BurnEditor {
  public BurnEditor(WindowRenderer.ManagerInterface manager,
                    IntPtr plugin,
                    Vessel vessel) {
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
                formatter        : value => FormatTimeSpan(
                                                TimeSpan.FromSeconds(value)));
    reference_frame_selector_ =
        new ReferenceFrameSelector(manager, plugin, ReferenceFrameChanged);
    plugin_ = plugin;
    vessel_ = vessel;
  }

  // Renders the |BurnEditor|.  Returns true if and only if the settings were
  // changed.
  public bool Render(bool enabled) {
    var old_skin = UnityEngine.GUI.skin;
    UnityEngine.GUI.skin = null;
    UnityEngine.GUILayout.BeginVertical();
    if (enabled) {
      reference_frame_selector_.RenderButton();
    } else {
      reference_frame_selector_.Hide();
    }
    bool changed = false;
    changed |= Δv_tangent_.Render(enabled);
    changed |= Δv_normal_.Render(enabled);
    changed |= Δv_binormal_.Render(enabled);
    changed |= initial_time_.Render(enabled);
    changed |= changed_reference_frame_;
    changed_reference_frame_ = false;
    UnityEngine.GUILayout.Label(
        text       : "(from end of previous burn)",
        options    : UnityEngine.GUILayout.Width(250));
    UnityEngine.GUILayout.EndVertical();
    UnityEngine.GUI.skin = old_skin;
    return changed && enabled;
  }

  public void Reset(Burn burn) {
    Δv_tangent_.value = burn.delta_v.x;
    Δv_normal_.value = burn.delta_v.y;
    Δv_binormal_.value = burn.delta_v.z;
    initial_time_.value = burn.initial_time;
    reference_frame_selector_.Reset(burn.frame);
  }

  public Burn Burn() {
    return new Burn{
        thrust_in_kilonewtons = thrust_in_kilonewtons_,
        specific_impulse_in_seconds_g0 = specific_impulse_in_seconds_g0_,
        frame = reference_frame_selector_.FrameParameters(),
        initial_time = initial_time_.value + flight_plan_initial_time,
        delta_v = new XYZ{x = Δv_tangent_.value,
                          y = Δv_normal_.value,
                          z = Δv_binormal_.value } };
  }

  public void ReferenceFrameChanged() {
    changed_reference_frame_ = true;
  }

  public void Close() {
    reference_frame_selector_.Dispose();
  }

  // Returns the equivalent of the .NET >= 4 format
  // span.ToString(@"ddd \d hh \h mm \m\i\n ss.FFF \s").
  private string FormatTimeSpan (TimeSpan span) {
     return span.Days.ToString("000") + " d " +
            span.Hours.ToString("00") + " h " +
            span.Minutes.ToString("00") + " min " +
            span.Seconds.ToString("00") + " s";
  }

  private DifferentialSlider Δv_tangent_;
  private DifferentialSlider Δv_normal_;
  private DifferentialSlider Δv_binormal_;
  private DifferentialSlider initial_time_;
  private ReferenceFrameSelector reference_frame_selector_;
  private double thrust_in_kilonewtons_;
  private double specific_impulse_in_seconds_g0_;

  private const double Log10ΔvLowerRate = -3.0;
  private const double Log10ΔvUpperRate = 3.5;
  private const double Log10TimeLowerRate = 0.0;
  private const double Log10TimeUpperRate = 7.0;

  // Not owned.
  private readonly IntPtr plugin_;
  private readonly Vessel vessel_;

  private bool changed_reference_frame_ = false;
}

}  // namespace ksp_plugin_adapter
}  // namespace principia
