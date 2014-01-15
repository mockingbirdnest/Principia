using KSP;
using NewtonianPhysics;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using UnityEngine;

namespace Principia {
  [KSPAddon(KSPAddon.Startup.Flight, false)]
  public class Core : MonoBehaviour {
    protected Rect windowPos;
    private Vector3d a;
    private Dictionary<string, Body> bodies;
    private bool fixOrbitalPositions = false;
    private Vector3d q;
    private NBodySystem system;
    private Vector3d v;

    private void drawGUI() {
      GUI.skin = HighLogic.Skin;
      windowPos = GUILayout.Window(1,
                                   windowPos,
                                   WindowGUI,
                                   "Traces of Various Descriptions",
                                   GUILayout.MinWidth(500));
    }
    private void FixedUpdate() {
      if (!FlightGlobals.ready) {
        return;
      }
      Vessel activeVessel = FlightGlobals.ActiveVessel;
      if (activeVessel != null) {
        a = new Vector3d {
          x = activeVessel.acceleration.x,
          y = activeVessel.acceleration.y,
          z = activeVessel.acceleration.z
        };
        q = new Vector3d {
          x = activeVessel.orbit.pos.x,
          y = activeVessel.orbit.pos.y,
          z = activeVessel.orbit.pos.z
        };
        v = new Vector3d {
          x = activeVessel.orbit.vel.x,
          y = activeVessel.orbit.vel.y,
          z = activeVessel.orbit.vel.z
        };
      }
      if (fixOrbitalPositions) {
        double UT = Planetarium.GetUniversalTime();
        QuaternionD rotation = Planetarium.Rotation;
#if TRACE
        print("Initiating simulation step...");
#endif
        system.Evolve(UT, 10);
        foreach (CelestialBody body in FlightGlobals.Bodies) {
          if (body.name != "Sun") {
#if TRACE
            print("Updating " + body.name + "...");
#endif
            Body secondary, primary;
            bodies.TryGetValue(body.name, out secondary);
            bodies.TryGetValue(body.referenceBody.name, out primary);
            Vector3d position = secondary.q.ToVector() - primary.q.ToVector();
            Vector3d velocity = secondary.v.ToVector() - primary.v.ToVector();
            Orbit original = body.orbitDriver.orbit;
            Orbit copy = new Orbit(original.inclination,
                                   original.eccentricity,
                                   original.semiMajorAxis,
                                   original.LAN,
                                   original.argumentOfPeriapsis,
                                   original.meanAnomalyAtEpoch,
                                   original.epoch,
                                   original.referenceBody);
            copy.UpdateFromStateVectors((rotation * position.xzy).xzy,
                                        (rotation * velocity.xzy).xzy,
                                        copy.referenceBody,
                                        UT);
            body.orbit.inclination = copy.inclination;
            body.orbit.eccentricity = copy.eccentricity;
            body.orbit.semiMajorAxis = copy.semiMajorAxis;
            body.orbit.LAN = copy.LAN;
            body.orbit.argumentOfPeriapsis = copy.argumentOfPeriapsis;
            body.orbit.meanAnomalyAtEpoch = copy.meanAnomalyAtEpoch;
            body.orbit.epoch = copy.epoch;
            body.orbit.referenceBody = copy.referenceBody;
            body.orbit.Init();
            body.orbit.UpdateFromUT(UT);
            body.CBUpdate();
          }
        }
        foreach (Vessel vessel in FlightGlobals.Vessels) {
          if (vessel.situation != Vessel.Situations.LANDED &&
            vessel.situation != Vessel.Situations.SPLASHED &&
            vessel.situation != Vessel.Situations.PRELAUNCH) {
#if TRACE
            print("Updating " + vessel.name + "...");
#endif
            Body secondary, primary;
            bodies.TryGetValue(vessel.id.ToString(), out secondary);
            bodies.TryGetValue(vessel.orbit.referenceBody.name, out primary);
            Vector3d position = secondary.q.ToVector() - primary.q.ToVector();
            Vector3d velocity = secondary.v.ToVector() - primary.v.ToVector();
            vessel.orbit.UpdateFromStateVectors((rotation * position.xzy).xzy,
                                                (rotation * velocity.xzy).xzy,
                                                vessel.orbit.referenceBody,
                                                UT);
          }
        }
      }
    }
    private void OnDestroy() {
      RenderingManager.RemoveFromPostDrawQueue(3, new Callback(drawGUI));
    }

    private void Start() {
      print("Principia: Start!");
      bodies = new Dictionary<string, Body>();
      RenderingManager.AddToPostDrawQueue(3, new Callback(drawGUI));
      if ((windowPos.x == 0) && (windowPos.y == 0)) {
        windowPos = new Rect(Screen.width / 2, Screen.height / 2, 10, 10);
      }
    }
    private void WindowGUI(int windowID) {
      GUIStyle mySty = new GUIStyle(GUI.skin.button);
      mySty.normal.textColor = mySty.focused.textColor = Color.white;
      mySty.hover.textColor = mySty.active.textColor = Color.yellow;
      mySty.onNormal.textColor
        = mySty.onFocused.textColor
        = mySty.onHover.textColor
        = mySty.onActive.textColor = Color.green;
      mySty.padding = new RectOffset(8, 8, 8, 8);

      GUILayout.BeginVertical();
      GUILayout.TextArea("q: " + q.x.ToString("F3") + ", "
                               + q.y.ToString("F3") + ", "
                               + q.z.ToString("F3") + " ("
                               + q.magnitude.ToString("F3") + ")",
                               GUILayout.ExpandWidth(true));
      GUILayout.TextArea("v: " + v.x.ToString("F3") + ", "
                               + v.y.ToString("F3") + ", "
                               + v.z.ToString("F3") + " ("
                               + v.magnitude.ToString("F3") + ")",
                               GUILayout.ExpandWidth(true));
      GUILayout.TextArea("a: " + a.x.ToString("F3") + ", "
                               + a.y.ToString("F3") + ", "
                               + a.z.ToString("F3") + " ("
                               + a.magnitude.ToString("F3") + ")",
                               GUILayout.ExpandWidth(true));
      GUILayout.TextArea("Planetarium rotation: "
                         + Planetarium.Rotation.w.ToString("F5") + " + "
                         + Planetarium.Rotation.x.ToString("F5") + " i + "
                         + Planetarium.Rotation.y.ToString("F5") + " j + "
                         + Planetarium.Rotation.z.ToString("F5") + " k",
                         GUILayout.ExpandWidth(true));
      if (GUILayout.Button(fixOrbitalPositions ? "Unlock" : "Lock",
                           mySty, GUILayout.ExpandWidth(true))) {
        fixOrbitalPositions = !fixOrbitalPositions;
        double UT = Planetarium.GetUniversalTime();
        bodies.Clear();
        foreach (CelestialBody body in FlightGlobals.Bodies) {
          print("Adding body " + body.name + "...");
          bodies.Add(body.name, body.ToBody());
        }
        foreach (Vessel vessel in FlightGlobals.Vessels) {
          if (vessel.situation == Vessel.Situations.ESCAPING ||
              vessel.situation == Vessel.Situations.ORBITING ||
              vessel.situation == Vessel.Situations.SUB_ORBITAL) {
            print("Adding vessel " + vessel.name + "...");
            bodies.Add(vessel.id.ToString(), vessel.ToBody());
          }
        }
        system = new NBodySystem(bodies.Values.ToArray(), UT);
      }
      if (fixOrbitalPositions) {
        GUILayout.TextArea("Simulating n-body physics...");
        Body kerbin;
        bodies.TryGetValue("Kerbin", out kerbin);
        GUILayout.TextArea("qKerbin: " + kerbin.q.x + ", "
                                       + kerbin.q.y + ", "
                                       + kerbin.q.z);
        GUILayout.TextArea("vKerbin: " + kerbin.v.x + ", "
                                       + kerbin.v.y + ", "
                                       + kerbin.v.z);
      }
      GUILayout.EndVertical();

      GUI.DragWindow(new Rect(0, 0, 10000, 20));
    }
  }
}