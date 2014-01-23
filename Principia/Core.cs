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
    private const double Day = 24 * Hour;
    private const double Hour = 60 * Minute;
    private const double Minute = 60 * Second;
    private const double PredictionLength = 1 * Day;
    private const double Second = 1;
    private const double Year = 365.25 * Day;
    private Vector3d a;
    private Dictionary<string, Body> bodies;
    private bool predictTrajectories = false;
    private Vector3d q;
    private bool simulate = false;
    private NBodySystem system;
    private Dictionary<string, LineRenderer> trajectories;
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
      double UT = Planetarium.GetUniversalTime();
      QuaternionD rotation = Planetarium.Rotation;
      if (simulate) {
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
      if (predictTrajectories) {
        system.AdvancePredictions(tmax: UT + 1 * Day,
                                  maxTimestep: 10 * Second,
                                  samplingPeriod: 100);
        if (MapView.MapIsEnabled) {
          PlanetariumCamera camera = (PlanetariumCamera)
            GameObject.FindObjectOfType(typeof(PlanetariumCamera));
          foreach (Vessel vessel in FlightGlobals.Vessels) {
            if (vessel.situation != Vessel.Situations.LANDED &&
              vessel.situation != Vessel.Situations.SPLASHED &&
              vessel.situation != Vessel.Situations.PRELAUNCH) {
              LineRenderer line;
              Body body;
              trajectories.TryGetValue(vessel.id.ToString(), out line);
              bodies.TryGetValue(vessel.id.ToString(), out body);
              line.SetVertexCount(body.predictedTrajectory.Count);
              line.enabled = true;
              for (int i = 0; i < body.predictedTrajectory.Count; ++i) {
                line.SetPosition(i, ScaledSpace.LocalToScaledSpace(
                  body.predictedTrajectory[i].q.ToVector().xzy
                  - body.q.ToVector().xzy + vessel.GetWorldPos3D()));
                line.SetWidth((0.01f * camera.Distance),
                              (0.01f * camera.Distance));
              }
            }
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
      trajectories = new Dictionary<string, LineRenderer>();
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
      if (GUILayout.Button(simulate ? "Switch to Keplerian" : "Switch to Newtonian",
                           mySty, GUILayout.ExpandWidth(true))) {
        simulate = !simulate;
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
      if (GUILayout.Button(
                    predictTrajectories ? "Stop Plotting" : "Plot Trajectories",
                    mySty, GUILayout.ExpandWidth(true))) {
        predictTrajectories = !predictTrajectories;
        double UT = Planetarium.GetUniversalTime();
        if (predictTrajectories) {
          trajectories.Clear();
          foreach (Vessel vessel in FlightGlobals.Vessels) {
            if (vessel.situation == Vessel.Situations.ESCAPING ||
                vessel.situation == Vessel.Situations.ORBITING ||
                vessel.situation == Vessel.Situations.SUB_ORBITAL) {
              print("Adding vessel " + vessel.name + " trajectory...");
              GameObject lineObject = new GameObject("Line");
              lineObject.layer = 31;
              LineRenderer line = lineObject.AddComponent<LineRenderer>();
              line.transform.parent = null;
              line.material = MapView.fetch.orbitLinesMaterial;
              line.SetColors(Color.blue, Color.red);
              line.useWorldSpace = true;
              trajectories.Add(vessel.id.ToString(), line);
            }
          }
          system.RecalculatePredictions(tmax: UT + 1 * Day,
                                        maxTimestep: 10 * Second,
                                        samplingPeriod: 100);
        }
      }
      if (simulate) {
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