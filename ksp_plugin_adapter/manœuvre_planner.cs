using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace principia {
namespace ksp_plugin_adapter {

internal class QuantityInputField {
  double value_;
  string text_;
  string old_text_;
  string label_;
  string status_;
  string unit_;

  public delegate string DoubleToString(double value);
  public delegate bool StringToDouble(string text,
                                      out double value,
                                      out string status);

  DoubleToString to_string_;
  StringToDouble to_double_;

  private static bool ParseDouble(string text,
                                  out double value,
                                  out string status) {
    if (double.TryParse(text, out value)) {
      status = "";
      return true;
    } else {
      status = "Invalid number";
      return false;
    }
  }

  public QuantityInputField(string label,
                            double value,
                            string unit) : this(label,
                            value,
                            unit,
                            (x => x.ToString()),
                            ParseDouble) {}


  public QuantityInputField(string label,
                            double value,
                            string unit,
                            DoubleToString to_string,
                            StringToDouble to_double) {
    to_string_ = to_string;
    to_double_ = to_double;
    label_ = label;
    value_ = value;
    text_ = to_string_(value);
    old_text_ = text_;
    unit_ = unit;
  }

  public void Render() {
    UnityEngine.GUILayout.Label(label_, UnityEngine.GUILayout.Width(75));
    text_ = UnityEngine.GUILayout.TextField(text_,
                                            UnityEngine.GUILayout.Width(200));
    UnityEngine.GUILayout.Label(unit_, UnityEngine.GUILayout.Width(75));
    UnityEngine.GUILayout.Label(status_, UnityEngine.GUILayout.Width(75));
    UnityEngine.GUILayout.EndHorizontal();
  }

  public double value {
    get {
      return value_;
    }
    set {
      value_ = value;
      text_ = to_string_(value_);
      old_text_ = text_;
      status_ = "";
    }
  }

  public void ReadValue() {
    if (old_text_ != text_) {
      double parsed_value;
      if (to_double_(text_, out parsed_value, out status_)) {
        value = parsed_value;
        text_ = to_string_(parsed_value);
        old_text_ = text_;
      } else {
        status_ = "";
      }
    }
  }
}

internal class ManœuvrePlanner {
  QuantityInputField thrust_input_;
  QuantityInputField initial_mass_input_;
  QuantityInputField specific_impulse_input_;
  QuantityInputField right_ascension_input_;
  QuantityInputField declination_input_;
  QuantityInputField Δv_input_;
  QuantityInputField initial_time_input_;

  double duration_ = double.NaN;
  double final_time_ = double.NaN;
  double final_mass_ = double.NaN;

  Vessel vessel_;
  string vessel_guid_;

  public ManœuvrePlanner previous { get; set; }

  internal ManœuvrePlanner(Vessel vessel, ManœuvrePlanner previous) {
    vessel_ = vessel;
    vessel_guid_ = vessel_.id.ToString();

    thrust_input_           = new QuantityInputField("F", 1, "N");
    initial_mass_input_     = new QuantityInputField("m_0", 100, "kg");
    specific_impulse_input_ = new QuantityInputField("I_sp", 1, "s g_0");
    right_ascension_input_  = new QuantityInputField("α", 0, "°");
    declination_input_      = new QuantityInputField("δ", 0, "°");
    Δv_input_               = new QuantityInputField("Δv", 1, "m/s");
    initial_time_input_     = new QuantityInputField("t_0", 0,
                                                     "s after launch",
                                                     InitialTimeToString,
                                                     ParseInitialTime);

    this.previous = previous;
  }

  internal void Load(IntPtr manœuvre) {
    thrust_input_.value = PrincipiaPluginAdapter.thrust(manœuvre);
    initial_mass_input_.value = PrincipiaPluginAdapter.initial_mass(manœuvre);
    specific_impulse_input_.value =
        PrincipiaPluginAdapter.specific_impulse_by_weight(manœuvre);
    right_ascension_input_.value =
        PrincipiaPluginAdapter.right_ascension(manœuvre);
    declination_input_.value = PrincipiaPluginAdapter.declination(manœuvre);
    Δv_input_.value = PrincipiaPluginAdapter.Δv(manœuvre);
    initial_time_input_.value = PrincipiaPluginAdapter.initial_time(manœuvre);

    duration_ = PrincipiaPluginAdapter.duration(manœuvre);
    final_time_ = PrincipiaPluginAdapter.duration(manœuvre);
  }

  internal void Render() {
    thrust_input_.Render();
    specific_impulse_input_.Render();
    if (UnityEngine.GUILayout.Button(
            text    : "Auto Engines",
            options : UnityEngine.GUILayout.Width(100)) &&
        FlightGlobals.fetch.activeVessel != null) {
      AutoEngines();
    }
    initial_mass_input_.Render();
    if (UnityEngine.GUILayout.Button(
            text    : "Auto Mass",
            options : UnityEngine.GUILayout.Width(100)) &&
        FlightGlobals.fetch.activeVessel != null) {
      if (previous != null) {
        initial_mass_input_.value = previous.final_mass_;
      } else {
        initial_mass_input_.value = vessel_.GetTotalMass() * 1000;
      }
    }
    right_ascension_input_.Render();
    declination_input_.Render();
    Δv_input_.Render();
    initial_time_input_.Render();
    UnityEngine.GUILayout.BeginHorizontal();
    var old_alignment = UnityEngine.GUI.skin.textArea.alignment;
    UnityEngine.GUI.skin.textArea.alignment =
        UnityEngine.TextAnchor.MiddleRight;
    UnityEngine.GUILayout.Label(text    : "Burn duration:",
                                options : UnityEngine.GUILayout.Width(150));
    UnityEngine.GUILayout.TextArea(
        text    : duration_.ToString() + " s",
        options : UnityEngine.GUILayout.Width(75));
    UnityEngine.GUILayout.EndHorizontal();
    UnityEngine.GUILayout.BeginHorizontal();
    UnityEngine.GUILayout.Label(text    : "Burn start countdown",
                                options : UnityEngine.GUILayout.Width(150));
    UnityEngine.GUILayout.TextArea(
        (initial_time_input_.value -
         Planetarium.GetUniversalTime()).ToString("+0.000;-0.000") + " s",
        UnityEngine.GUILayout.Width(75));
    UnityEngine.GUILayout.Label(text    : "Burn end countdown",
                                options : UnityEngine.GUILayout.Width(150));
    UnityEngine.GUILayout.TextArea(
        (final_time_ -
         Planetarium.GetUniversalTime()).ToString("+0.000;-0.000") + " s",
        UnityEngine.GUILayout.Width(75));
    UnityEngine.GUILayout.EndHorizontal();
    UnityEngine.GUI.skin.textArea.alignment = old_alignment;
  }

  internal IntPtr ComputeManœuvre() {
    thrust_input_.ReadValue();
    specific_impulse_input_.ReadValue();
    initial_mass_input_.ReadValue();
    right_ascension_input_.ReadValue();
    declination_input_.ReadValue();
    Δv_input_.ReadValue();
    initial_time_input_.ReadValue();

    IntPtr manœuvre = PrincipiaPluginAdapter.NewManœuvreIspByWeight(
        thrust_input_.value,
        initial_time_input_.value,
        specific_impulse_input_.value,
        right_ascension_input_.value,
        declination_input_.value);
    PrincipiaPluginAdapter.set_Δv(manœuvre, Δv_input_.value);
    PrincipiaPluginAdapter.set_initial_time(manœuvre, initial_time_input_.value);

    Load(manœuvre);
    return manœuvre;
  }

  private void AutoEngines() {
    ModuleEngines[] active_engines =
        (from part in vessel_.parts
         select (from PartModule module in part.Modules
                 where module is ModuleEngines &&
                       (module as ModuleEngines).EngineIgnited
                 select module as ModuleEngines))
            .SelectMany(x => x)
            .ToArray();
    Vector3d reference_direction = vessel_.ReferenceTransform.up;
    double[] thrusts =
        (from engine in active_engines
         select engine.maxThrust * 1000 *
                (from transform in engine.thrustTransforms
                 select Vector3d.Dot(reference_direction,
                                     -transform.forward)).Average())
            .ToArray();
    double total_thrust = thrusts.Sum();
    thrust_input_.value = total_thrust;

    // This would use zip if we had 4.0 or later.  We loop for now.
    double Σ_f_over_i_sp = 0;
    for (int i = 0; i < active_engines.Count(); ++i) {
      Σ_f_over_i_sp += thrusts[i] /
                       active_engines[i].atmosphereCurve.Evaluate(0);
    }
    specific_impulse_input_.value = total_thrust / Σ_f_over_i_sp;
  }

  private bool ParseInitialTime(string text,
                                out double value,
                                out string status) {
    if (double.TryParse(text, out value)) {
      status = "";
      value += vessel_.launchTime;
      if (value <= Planetarium.GetUniversalTime()) {
        status = "Time is in the past";
        return false;
      } else if (previous != null && value <= previous.final_time_) {
        status = "Time before end of previous manœuvre";
        return false;
      } else {
        status = "";
        return true;
      }
    } else {
      status = "Invalid number";
      return false;
    }
  }

  private string InitialTimeToString(double value) {
    return (value - vessel_.launchTime).ToString();
  }
}

}  // namespace ksp_plugin_adapter
}  // namespace principia
