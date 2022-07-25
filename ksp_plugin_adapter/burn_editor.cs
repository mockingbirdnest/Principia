using System;
using System.Linq;
using System.Text.RegularExpressions;
using KSP.Localization;

namespace principia {
namespace ksp_plugin_adapter {

class BurnEditor : ScalingRenderer {
  public BurnEditor(PrincipiaPluginAdapter adapter,
                    Vessel vessel,
                    double initial_time,
                    int index,
                    Func<int, BurnEditor> get_burn_at_index) {
    adapter_ = adapter;
    vessel_ = vessel;
    initial_time_ = initial_time;
    this.index = index;
    get_burn_at_index_ = get_burn_at_index;
    Δv_tangent_ = new DifferentialSlider(
        label            : L10N.CacheFormat("#Principia_BurnEditor_ΔvTangent"),
        unit             : L10N.CacheFormat("#Principia_BurnEditor_SpeedUnit"),
        min_value        : -max_Δv_component,
        max_value        : max_Δv_component,
        log10_lower_rate : log10_Δv_lower_rate,
        log10_upper_rate : log10_Δv_upper_rate,
        formatter        : FormatΔvComponent,
        alignment        : UnityEngine.TextAnchor.MiddleLeft,
        text_colour      : Style.Tangent);
    Δv_normal_ = new DifferentialSlider(
        label            : L10N.CacheFormat("#Principia_BurnEditor_ΔvNormal"),
        unit             : L10N.CacheFormat("#Principia_BurnEditor_SpeedUnit"),
        min_value        : -max_Δv_component,
        max_value        : max_Δv_component,
        log10_lower_rate : log10_Δv_lower_rate,
        log10_upper_rate : log10_Δv_upper_rate,
        formatter        : FormatΔvComponent,
        alignment        : UnityEngine.TextAnchor.MiddleLeft,
        text_colour      : Style.Normal);
    Δv_binormal_ = new DifferentialSlider(
        label            : L10N.CacheFormat("#Principia_BurnEditor_ΔvBinormal"),
        unit             : L10N.CacheFormat("#Principia_BurnEditor_SpeedUnit"),
        min_value        : -max_Δv_component,
        max_value        : max_Δv_component,
        log10_lower_rate : log10_Δv_lower_rate,
        log10_upper_rate : log10_Δv_upper_rate,
        formatter        : FormatΔvComponent,
        alignment        : UnityEngine.TextAnchor.MiddleLeft,
        text_colour      : Style.Binormal);
    previous_coast_duration_ = new DifferentialSlider(
        label            : L10N.CacheFormat("#Principia_BurnEditor_InitialTime"),
        unit             : null,
        log10_lower_rate : log10_time_lower_rate,
        log10_upper_rate : log10_time_upper_rate,
        // We cannot have a coast of length 0, so let's make it very
        // short: that will be indistinguishable.
        zero_value       : 0.001,
        min_value        : 0,
        formatter        : FormatPreviousCoastDuration,
        parser           : TryParsePreviousCoastDuration,
        field_width      : 7){
        value            = initial_time_ - time_base
    };
    reference_frame_selector_ = new ReferenceFrameSelector(
        adapter_,
        ReferenceFrameChanged,
        L10N.CacheFormat("#Principia_BurnEditor_ManœuvringFrame"));
    reference_frame_selector_.SetFrameParameters(
        adapter_.plotting_frame_selector_.FrameParameters());
    ComputeEngineCharacteristics();
  }

  public enum Event {
    None,
    Changed,
    Deleted,
    Minimized,
    Maximized,
  }

  // Renders the |BurnEditor|.  Returns true if and only if the settings were
  // changed.
  public Event Render(
      string header,
      bool anomalous,
      double burn_final_time,
      double? orbital_period) {
    bool changed = false;
    previous_coast_duration_.max_value = burn_final_time - time_base;
    using (new UnityEngine.GUILayout.HorizontalScope()) {
      if (UnityEngine.GUILayout.Button(
              minimized ? "+" : "−", GUILayoutWidth(1))) {
        minimized = !minimized;
        return minimized ? Event.Minimized : Event.Maximized;
      }
      UnityEngine.GUILayout.Label(
          minimized ? L10N.CacheFormat("#Principia_BurnEditor_MinimizedHeader",
                                       header,
                                       Δv().ToString("0.000"))
                    : header);
      string info = "";
      if (!minimized &&
          !reference_frame_selector_.FrameParameters().Equals(
              adapter_.plotting_frame_selector_.FrameParameters())) {
        info = L10N.CacheFormat(
            "#Principia_BurnEditor_Info_InconsistentFrames");
      }
      UnityEngine.GUILayout.Label(info, Style.Info(UnityEngine.GUI.skin.label));
      if (UnityEngine.GUILayout.Button(
              L10N.CacheFormat("#Principia_BurnEditor_Delete"),
              GUILayoutWidth(2))) {
        return Event.Deleted;
      }
    }
    if (minimized) {
      return Event.None;
    }
    using (new UnityEngine.GUILayout.VerticalScope()) {
      // When we are first rendered, the |initial_mass_in_tonnes_| will just have
      // been set.  If we have fallen back to instant impulse, we should use this
      // mass to set the thrust.
      if (first_time_rendering_) {
        first_time_rendering_ = false;
        changed = true;
        engine_warning_ = "";
        ComputeEngineCharacteristics();
      }

      // The frame selector is disabled for an anomalous manœuvre as is has no
      // effect.
      if (anomalous) {
        reference_frame_selector_.Hide();
      } else {
        using (new UnityEngine.GUILayout.HorizontalScope()) {
          if (UnityEngine.GUILayout.Button(
              L10N.CacheFormat("#Principia_BurnEditor_ActiveEngines"))) {
            engine_warning_ = "";
            ComputeEngineCharacteristics();
            changed = true;
          } else if (UnityEngine.GUILayout.Button(
              L10N.CacheFormat("#Principia_BurnEditor_ActiveRCS"))) {
            engine_warning_ = "";
            ComputeRCSCharacteristics();
            changed = true;
          } else if (UnityEngine.GUILayout.Button(
              L10N.CacheFormat("#Principia_BurnEditor_InstantImpulse"))) {
            engine_warning_ = "";
            UseTheForceLuke();
            changed = true;
          }
        }
        reference_frame_selector_.RenderButton();
      }
      if (is_inertially_fixed_ !=
          UnityEngine.GUILayout.Toggle(
              is_inertially_fixed_,
              L10N.CacheFormat("#Principia_BurnEditor_InertiallyFixed"))) {
        changed = true;
        is_inertially_fixed_ = !is_inertially_fixed_;
      }
      changed |= changed_reference_frame_;

      // The Δv controls are disabled for an anomalous manœuvre as they have no
      // effect.
      changed |= Δv_tangent_.Render(enabled : !anomalous);
      changed |= Δv_normal_.Render(enabled : !anomalous);
      changed |= Δv_binormal_.Render(enabled : !anomalous);
      {
        var render_time_base = time_base;
        previous_coast_duration_.value_if_different =
            initial_time_ - render_time_base;
        // The duration of the previous coast is always enabled as it can make
        // a manœuvre non-anomalous.
        if (previous_coast_duration_.Render(enabled : true)) {
          changed = true;
          initial_time_ = previous_coast_duration_.value + render_time_base;
        }
      }
      using (new UnityEngine.GUILayout.HorizontalScope()) {
        UnityEngine.GUILayout.Label("", GUILayoutWidth(3));
        if (decrement_revolution == null) {
          PrincipiaPluginAdapter.LoadTextureOrDie(out decrement_revolution,
                                                  "decrement_revolution.png");
        }
        if (increment_revolution == null) {
          PrincipiaPluginAdapter.LoadTextureOrDie(out increment_revolution,
                                                  "increment_revolution.png");
        }
        if (orbital_period is double period) {
          if (UnityEngine.GUILayout.Button(
                  new UnityEngine.GUIContent(
                      decrement_revolution,
                      L10N.CacheFormat(
                          "#Principia_BurnEditor_DecrementRevolution")),
                  GUILayoutWidth(1))) {
            changed = true;
            initial_time_ -= period;
          }
          UnityEngine.GUILayout.Space(Width(5));
          if (UnityEngine.GUILayout.Button(
                  new UnityEngine.GUIContent(
                      increment_revolution,
                      L10N.CacheFormat(
                          "#Principia_BurnEditor_IncrementRevolution")),
                  GUILayoutWidth(1))) {
            changed = true;
            initial_time_ += period;
          }
        } else {
          UnityEngine.GUILayout.Button("", GUILayoutWidth(1));
          UnityEngine.GUILayout.Space(Width(5));
          UnityEngine.GUILayout.Button("", GUILayoutWidth(1));
        }
        UnityEngine.GUILayout.Label(
            index == 0
                ? L10N.CacheFormat(
                    "#Principia_BurnEditor_TimeBase_StartOfFlightPlan")
                : L10N.CacheFormat("#Principia_BurnEditor_TimeBase_EndOfManœuvre",
                                   index),
            style : new UnityEngine.GUIStyle(UnityEngine.GUI.skin.label){
                alignment = UnityEngine.TextAnchor.UpperLeft
            });
      }
      using (new UnityEngine.GUILayout.HorizontalScope()) {
        UnityEngine.GUILayout.Label(
            L10N.CacheFormat("#Principia_BurnEditor_Δv",
                             Δv().ToString("0.000")),
            GUILayoutWidth(8));
        UnityEngine.GUILayout.Label(L10N.CacheFormat(
                                        "#Principia_BurnEditor_Duration",
                                        duration_.ToString("0.0")));
      }
      UnityEngine.GUILayout.Label(engine_warning_,
                                  Style.Warning(UnityEngine.GUI.skin.label));
      changed_reference_frame_ = false;
    }
    return changed ? Event.Changed : Event.None;
  }

  public double Δv() {
    return new Vector3d{
        x = Δv_tangent_.value,
        y = Δv_normal_.value,
        z = Δv_binormal_.value
    }.magnitude;
  }

  public void Reset(NavigationManoeuvre manœuvre) {
    Burn burn = manœuvre.burn;
    Δv_tangent_.value = burn.delta_v.x;
    Δv_normal_.value = burn.delta_v.y;
    Δv_binormal_.value = burn.delta_v.z;
    initial_time_ = burn.initial_time;
    reference_frame_selector_.SetFrameParameters(burn.frame);
    is_inertially_fixed_ = burn.is_inertially_fixed;
    duration_ = manœuvre.duration;
    initial_mass_in_tonnes_ = manœuvre.initial_mass_in_tonnes;
  }

  public Burn Burn() {
    return new Burn{
        thrust_in_kilonewtons = thrust_in_kilonewtons_,
        specific_impulse_in_seconds_g0 = specific_impulse_in_seconds_g0_,
        frame = reference_frame_selector_.FrameParameters(),
        initial_time = initial_time_,
        delta_v = new XYZ{
            x = Δv_tangent_.value,
            y = Δv_normal_.value,
            z = Δv_binormal_.value
        },
        is_inertially_fixed = is_inertially_fixed_
    };
  }

  public void ReferenceFrameChanged(NavigationFrameParameters? parameters,
                                    Vessel target_vessel) {
    changed_reference_frame_ = true;
  }

  public void Close() {
    reference_frame_selector_.DisposeWindow();
  }

  private void ComputeEngineCharacteristics() {
    ModuleEngines[] active_engines =
        (from part in vessel_.parts
         select (from PartModule module in part.Modules
                 where (module as ModuleEngines)?.EngineIgnited == true
                 select module as ModuleEngines)).SelectMany(x => x).ToArray();
    Vector3d reference_direction = vessel_.ReferenceTransform.up;
    double[] thrusts =
        (from engine in active_engines
         select
             engine.MaxThrustOutputVac(useThrustLimiter: true) *
             (from transform in engine.thrustTransforms
              select Math.Max(0,
                              Vector3d.Dot(reference_direction,
                                           -transform.forward))).Average()).
        ToArray();
    thrust_in_kilonewtons_ = thrusts.Sum();

    // This would use zip if we had 4.0 or later.  We loop for now.
    double Σ_f_over_i_sp = 0;
    for (int i = 0; i < active_engines.Length; ++i) {
      Σ_f_over_i_sp +=
          thrusts[i] / active_engines[i].atmosphereCurve.Evaluate(0);
    }
    specific_impulse_in_seconds_g0_ = thrust_in_kilonewtons_ / Σ_f_over_i_sp;

    // If there are no engines, fall back onto RCS.
    if (thrust_in_kilonewtons_ == 0) {
      engine_warning_ +=
          L10N.CacheFormat("#Principia_BurnEditor_Warning_NoActiveEngines");
      ComputeRCSCharacteristics();
    }
  }

  private void ComputeRCSCharacteristics() {
    ModuleRCS[] active_rcs = (from part in vessel_.parts
                              select (from PartModule module in part.Modules
                                      where module is ModuleRCS module_rcs &&
                                            module_rcs.rcsEnabled
                                      select module as ModuleRCS)).
        SelectMany(x => x).ToArray();
    Vector3d reference_direction = vessel_.ReferenceTransform.up;
    // NOTE(egg): NathanKell informs me that in >= 1.0.5, RCS has a useZaxis
    // property, that controls whether they thrust -up or -forward.  The madness
    // keeps piling up.
    double[] thrusts = (from engine in active_rcs
                        select engine.thrusterPower *
                               (from transform in engine.thrusterTransforms
                                select Math.Max(0,
                                                Vector3d.Dot(
                                                    reference_direction,
                                                    -transform.up))).Sum()).
        ToArray();
    thrust_in_kilonewtons_ = thrusts.Sum();

    // This would use zip if we had 4.0 or later.  We loop for now.
    double Σ_f_over_i_sp = 0;
    for (int i = 0; i < active_rcs.Length; ++i) {
      Σ_f_over_i_sp +=
          thrusts[i] / active_rcs[i].atmosphereCurve.Evaluate(0);
    }
    specific_impulse_in_seconds_g0_ = thrust_in_kilonewtons_ / Σ_f_over_i_sp;

    // If RCS provides no thrust, model a virtually instant burn.
    if (thrust_in_kilonewtons_ == 0) {
      engine_warning_ +=
          L10N.CacheFormat("#Principia_BurnEditor_Warning_NoActiveRCS");
      UseTheForceLuke();
    }
  }

  private string FormatΔvComponent(double metres_per_second) {
    // The granularity of Instant in 1950.
    double dt = 2.3841857910156250e-7;  // 2⁻²² s.
    double initial_acceleration =
        thrust_in_kilonewtons_ / initial_mass_in_tonnes_;
    double Isp = specific_impulse_in_seconds_g0_ * 9.80665;
    // Allow an increment in speed which, if applied to the norm of the
    // velocity, will always result in a change in the duration of the
    // integrated burn.
    double dv = dt * Math.Exp(Δv() / Isp) * initial_acceleration;
    // Allow no less than 1 nm/s to ensure that the number fits in the field,
    // and that the increment remains representable at any reasonable speed in
    // the solar system (that second criterion would be satisfied at 0.1 nm/s,
    // but the first would not).
    int fractional_digits =
        Math.Max(0, Math.Min((int)Math.Floor(-Math.Log10(dv)), 9));
    string unsigned_format = "00,000." + new string('0', fractional_digits);
    return Regex.Replace(
        Regex.Replace(
        metres_per_second.ToString($"+{unsigned_format};−{unsigned_format}",
                                   Culture.culture),
        "^[+−][0']{1,5}",
        match => match.Value.Replace('0', figure_space)),
        // Add grouping marks to the fractional part.
        @"\d{3}(?=\d)",
        match => match.Value + "'");
  }

  private void UseTheForceLuke() {
    // The burn can last at most (9.80665 / scale) s.
    const double scale = 1;
    // This, together with |scale = 1|, ensures that, when |initial_time| is
    // less than 2 ** 32 s, |Δv(initial_time + duration)| does not overflow if
    // Δv is less than 100 km/s, and that |initial_time + duration| does not
    // fully cancel if Δv is more than 1 mm/s.
    // TODO(egg): Before the C* release, add a persisted flag to indicate to the
    // user that we are not using the craft's engines (we can also use that
    // flag to remember whether the burn was created for active engines or
    // active RCS).
    const double range = 1000;
    thrust_in_kilonewtons_ = initial_mass_in_tonnes_ * range * scale;
    specific_impulse_in_seconds_g0_ = range;
  }

  internal string FormatPreviousCoastDuration(double seconds) {
    return new PrincipiaTimeSpan(seconds).FormatPositive(
        with_leading_zeroes: true,
        with_seconds: true,
        iau_style: true,
        fractional_second_digits: 6);
  }

  internal bool TryParsePreviousCoastDuration(string text, out double value) {
    value = 0;
    if (!PrincipiaTimeSpan.TryParse(text,
                                    out PrincipiaTimeSpan ts)) {
      return false;
    }
    value = ts.total_seconds;
    return true;
  }

  private double time_base => previous_burn?.final_time ??
                              plugin.FlightPlanGetInitialTime(
                                  vessel_.id.ToString());

  public double initial_time => initial_time_;
  public double final_time => initial_time_ + duration_;

  private IntPtr plugin => adapter_.Plugin();
  public int index { private get; set; }
  public bool minimized { private get; set; } = true;
  private BurnEditor previous_burn => get_burn_at_index_(index - 1);

  private bool is_inertially_fixed_;
  private readonly DifferentialSlider Δv_tangent_;
  private readonly DifferentialSlider Δv_normal_;
  private readonly DifferentialSlider Δv_binormal_;
  private readonly DifferentialSlider previous_coast_duration_;
  private readonly ReferenceFrameSelector reference_frame_selector_;
  private double thrust_in_kilonewtons_;
  private double specific_impulse_in_seconds_g0_;
  private double duration_;
  private double initial_mass_in_tonnes_;
  private double initial_time_;

  private bool first_time_rendering_ = true;

  private const double log10_Δv_lower_rate = -3.0;
  private const double log10_Δv_upper_rate = 3.5;
  private const double log10_time_lower_rate = 0.0;
  private const double log10_time_upper_rate = 7.0;
  private const double max_Δv_component = 99_999.999_999_999;

  // Not owned.
  private readonly Vessel vessel_;
  private readonly Func<int, BurnEditor> get_burn_at_index_;
  private readonly PrincipiaPluginAdapter adapter_;

  private bool changed_reference_frame_ = false;
  private string engine_warning_ = "";
  
  private static UnityEngine.Texture decrement_revolution;
  private static UnityEngine.Texture increment_revolution;
  private const char figure_space = '\u2007';
}

}  // namespace ksp_plugin_adapter
}  // namespace principia
