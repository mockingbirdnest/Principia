using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace principia {
namespace ksp_plugin_adapter {

class ManœuvrePlanner {
  string thrust_text_ = "1.0";
  string initial_mass_text_= "100.0";
  string specific_impulse_text_ = "1.0";
  string right_ascension_text_ = "0.0";
  string declination_text_ = "0.0";
  string Δv_text_ = "1.0";
  string initial_time_text_ = "10.0";

  double thrust_;
  double initial_mass_;
  double specific_impulse_;
  double right_ascension_;
  double declination_;
  double duration_;
  double Δv_;
  double initial_time_;

  internal void ReadManœuvre(IntPtr manœuvre) {
    thrust_ = PrincipiaPluginAdapter.thrust(manœuvre);
    initial_mass_ = PrincipiaPluginAdapter.initial_mass(manœuvre);
    specific_impulse_ =
        PrincipiaPluginAdapter.specific_impulse_by_weight(manœuvre);
    right_ascension_ = PrincipiaPluginAdapter.right_ascension(manœuvre);
    declination_ = PrincipiaPluginAdapter.declination(manœuvre);
    duration_ = PrincipiaPluginAdapter.duration(manœuvre);
    Δv_ = PrincipiaPluginAdapter.Δv(manœuvre);
    initial_time_ = PrincipiaPluginAdapter.initial_time(manœuvre);
    thrust_text_ = thrust_.ToString();
    initial_mass_text_ = initial_mass_.ToString();
    specific_impulse_text_ = specific_impulse_.ToString();
    declination_text_ = declination_.ToString();
    Δv_text_ = declination_.ToString();
    initial_time_text_
  }

  internal void Render() {
    UnityEngine.GUILayout.Label(text    : "F = ",
                                options : UnityEngine.GUILayout.Width(75));
    thrust_text_ =
        UnityEngine.GUILayout.TextField(thrust_text_,
                                        UnityEngine.GUILayout.Width(75));
    UnityEngine.GUILayout.Label(text : "N");
    UnityEngine.GUILayout.EndHorizontal();
    UnityEngine.GUILayout.BeginHorizontal();
    UnityEngine.GUILayout.Label(text    : "m_0 = ",
                                options : UnityEngine.GUILayout.Width(75));
    initial_mass_text_ =
        UnityEngine.GUILayout.TextField(initial_mass_text_,
                                        UnityEngine.GUILayout.Width(75));
    UnityEngine.GUILayout.Label(text : "kg");
    UnityEngine.GUILayout.EndHorizontal();
    UnityEngine.GUILayout.BeginHorizontal();
    UnityEngine.GUILayout.Label(text    : "I_sp = ",
                                options : UnityEngine.GUILayout.Width(75));
    specific_impulse_text_ =
        UnityEngine.GUILayout.TextField(specific_impulse_text_,
                                        UnityEngine.GUILayout.Width(75));
    UnityEngine.GUILayout.Label(text : "s g_0");
    UnityEngine.GUILayout.EndHorizontal();
    UnityEngine.GUILayout.BeginHorizontal();
    UnityEngine.GUILayout.Label(text    : "α = ",
                                options : UnityEngine.GUILayout.Width(75));
    right_ascension_text_ =
        UnityEngine.GUILayout.TextField(right_ascension_text_,
                                        UnityEngine.GUILayout.Width(75));
    UnityEngine.GUILayout.Label(text : "°");
    UnityEngine.GUILayout.EndHorizontal();
    UnityEngine.GUILayout.BeginHorizontal();
    UnityEngine.GUILayout.Label(text    : "δ = ",
                                options : UnityEngine.GUILayout.Width(75));
    declination_text_ =
        UnityEngine.GUILayout.TextField(declination_text_,
                                        UnityEngine.GUILayout.Width(75));
    UnityEngine.GUILayout.Label(text : "°");
    UnityEngine.GUILayout.EndHorizontal();
    UnityEngine.GUILayout.BeginHorizontal();
    UnityEngine.GUILayout.Label(text    : "Δv = ",
                                options : UnityEngine.GUILayout.Width(75));
    Δv_text_ = UnityEngine.GUILayout.TextField(Δv_text_,
                                         UnityEngine.GUILayout.Width(75));
    UnityEngine.GUILayout.Label(text : "m/s");
    UnityEngine.GUILayout.EndHorizontal();
    UnityEngine.GUILayout.BeginHorizontal();
    UnityEngine.GUILayout.Label(text    : "t_0 = ",
                                options : UnityEngine.GUILayout.Width(75));
    initial_time_text_ =
        UnityEngine.GUILayout.TextField(initial_time_text_,
                                        UnityEngine.GUILayout.Width(75));
    UnityEngine.GUILayout.Label(text : "s after launch.");
    UnityEngine.GUILayout.EndHorizontal();
  }
}

}  // namespace ksp_plugin_adapter
}  // namespace principia
