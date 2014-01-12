using KSP;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using UnityEngine;

namespace Principia {

  [KSPAddon(KSPAddon.Startup.Flight, false)]
  public class Core : MonoBehaviour {
    private Vector3d v;
    private Vector3d q;
    private Vector3d a;
    protected Rect windowPos;
    private bool fixOrbitalPositions = false;
    private Dictionary<string, Vector3d> positions;
    private Dictionary<string, Vector3d> velocities;

    private void WindowGUI(int windowID) {
      GUIStyle mySty = new GUIStyle(GUI.skin.button);
      mySty.normal.textColor = mySty.focused.textColor = Color.white;
      mySty.hover.textColor = mySty.active.textColor = Color.yellow;
      mySty.onNormal.textColor = mySty.onFocused.textColor = mySty.onHover.textColor = mySty.onActive.textColor = Color.green;
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
      if (GUILayout.Button(fixOrbitalPositions ? "Unlock" : "Lock", mySty, GUILayout.ExpandWidth(true))) {
        fixOrbitalPositions = !fixOrbitalPositions;
        double UT = Planetarium.GetUniversalTime();
        positions.Clear();
        velocities.Clear();
        foreach (CelestialBody body in FlightGlobals.Bodies) {
          if (body.name == "Sun") {
            continue;
          }
          positions.Add(body.name, (QuaternionD.Inverse(Planetarium.Rotation) * body.orbit.getRelativePositionAtUT(UT).xzy).xzy);
          velocities.Add(body.name, (QuaternionD.Inverse(Planetarium.Rotation) * body.orbit.getOrbitalVelocityAtUT(UT).xzy).xzy);
        }
        foreach (Vessel vessel in FlightGlobals.Vessels) {
          if (vessel.situation == Vessel.Situations.ESCAPING ||
              vessel.situation == Vessel.Situations.ORBITING ||
              vessel.situation == Vessel.Situations.SUB_ORBITAL) {
            positions.Add(vessel.id.ToString(), (QuaternionD.Inverse(Planetarium.Rotation) * vessel.orbit.getRelativePositionAtUT(UT).xzy).xzy);
            velocities.Add(vessel.id.ToString(), (QuaternionD.Inverse(Planetarium.Rotation) * vessel.orbit.getOrbitalVelocityAtUT(UT).xzy).xzy);
          }
        }
      }
      if (fixOrbitalPositions) {
        GUILayout.TextArea("Fixing orbital positions...");
        Vector3d position;
        positions.TryGetValue("Mun", out position);
        Vector3d velocity;
        velocities.TryGetValue("Mun", out velocity);
        GUILayout.TextArea("qMun: " + position.x + ", " + position.y + ", " + position.z);
        GUILayout.TextArea("vMun: " + velocity.x + ", " + velocity.y + ", " + velocity.z);
      }
      GUILayout.EndVertical();

      //DragWindow makes the window draggable. The Rect specifies which part of the window it can by dragged by, and is
      //clipped to the actual boundary of the window. You can also pass no argument at all and then the window can by
      //dragged by any part of it. Make sure the DragWindow command is AFTER all your other GUI input stuff, or else
      //it may "cover up" your controls and make them stop responding to the mouse.
      GUI.DragWindow(new Rect(0, 0, 10000, 20));
    }

    private void drawGUI() {
      GUI.skin = HighLogic.Skin;
      windowPos = GUILayout.Window(1, windowPos, WindowGUI, "Traces", GUILayout.MinWidth(500));
    }

    private void Start() {
      print("Miscellaneous Traces: Start!");
      positions = new Dictionary<string, Vector3d>();
      velocities = new Dictionary<string, Vector3d>();
      RenderingManager.AddToPostDrawQueue(3, new Callback(drawGUI));
      if ((windowPos.x == 0) && (windowPos.y == 0))//windowPos is used to position the GUI window, lets set it in the center of the screen
			{
        windowPos = new Rect(Screen.width / 2, Screen.height / 2, 10, 10);
      }
    }

    private void OnDestroy() {
      RenderingManager.RemoveFromPostDrawQueue(3, new Callback(drawGUI));
    }

    private void FixedUpdate() {
      if (!FlightGlobals.ready) {
        return;
      }
      Vessel activeVessel = FlightGlobals.ActiveVessel;
      if (activeVessel != null) {
        a = new Vector3d(activeVessel.acceleration.x, activeVessel.acceleration.y, activeVessel.acceleration.z);
        q = new Vector3d(activeVessel.orbit.pos.x, activeVessel.orbit.pos.y, activeVessel.orbit.pos.z);
        v = new Vector3d(activeVessel.orbit.vel.x, activeVessel.orbit.vel.y, activeVessel.orbit.vel.z);
      }
      if (fixOrbitalPositions) {
        foreach (CelestialBody body in FlightGlobals.Bodies) {
          if (body.name == "Sun") {
            continue;
          }
          Vector3d position;
          positions.TryGetValue(body.name, out position);
          Vector3d velocity;
          velocities.TryGetValue(body.name, out velocity);
          Orbit original = body.orbitDriver.orbit;
          Orbit copy = new Orbit(original.inclination,
                                 original.eccentricity,
                                 original.semiMajorAxis,
                                 original.LAN,
                                 original.argumentOfPeriapsis,
                                 original.meanAnomalyAtEpoch,
                                 original.epoch,
                                 original.referenceBody);
          copy.UpdateFromStateVectors((Planetarium.Rotation * position.xzy).xzy,
                                      (Planetarium.Rotation * velocity.xzy).xzy,
                                      copy.referenceBody,
                                      Planetarium.GetUniversalTime());
          body.orbit.inclination = copy.inclination;
          body.orbit.eccentricity = copy.eccentricity;
          body.orbit.semiMajorAxis = copy.semiMajorAxis;
          body.orbit.LAN = copy.LAN;
          body.orbit.argumentOfPeriapsis = copy.argumentOfPeriapsis;
          body.orbit.meanAnomalyAtEpoch = copy.meanAnomalyAtEpoch;
          body.orbit.epoch = copy.epoch;
          body.orbit.referenceBody = copy.referenceBody;
          body.orbit.Init();
          body.orbit.UpdateFromUT(Planetarium.GetUniversalTime());
          body.CBUpdate();
          /*   body.orbit.UpdateFromStateVectors((Planetarium.Rotation * position.xzy).xzy,
                                         (Planetarium.Rotation * velocity.xzy).xzy,
                                         body.referenceBody,
                                         Planetarium.GetUniversalTime());*/
        }
        foreach (Vessel vessel in FlightGlobals.Vessels) {
          if (vessel.situation != Vessel.Situations.LANDED &&
            vessel.situation != Vessel.Situations.SPLASHED &&
            vessel.situation != Vessel.Situations.PRELAUNCH) {
            Vector3d position;
            positions.TryGetValue(vessel.id.ToString(), out position);
            Vector3d velocity;
            velocities.TryGetValue(vessel.id.ToString(), out velocity);
            vessel.orbit.UpdateFromStateVectors((Planetarium.Rotation * position.xzy).xzy,
                                                (Planetarium.Rotation * velocity.xzy).xzy,
                                                vessel.orbit.referenceBody,
                                                Planetarium.GetUniversalTime());
          }
        }
      }
    }
  }
}