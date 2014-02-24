using KSP;
using NewtonianPhysics;
using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using UnityEngine;

namespace Principia {
  [KSPAddon(KSPAddon.Startup.Flight, false)]
  public class Core : MonoBehaviour {
    // TODO(robin): This class is a mess. The calculations, both for renderering
    // and moving vessels around, should be done elsewhere.
    protected Rect mainWindowPosition;
    protected Rect referenceFrameWindowPosition;
    private const double Day = 24 * Hour;
    private const double Hour = 60 * Minute;
    private const double Minute = 60 * Second;
    private const double PredictionLength = 1 * Day;
    private const double Second = 1;
    private const double Year = 365.25 * Day;
    private Vector3d activeVesselProperAcceleration;
    private Vector3d activeVesselVelocity;
    private Dictionary<string, Body> bodies;
    private Vector3d geometricAcceleration;
    private bool lastFrameWasPhysical = false;
    private MapRenderer mapRenderer;
    private bool predictTrajectories = false;
    private LineRenderer properAccelerationVisualisation;
    private CelestialBody referenceBody;
    private ReferenceFrameType referenceFrameType;
    private bool simulate = false;
    private Vector3d[] standardBasis = { new Vector3d { x = 1, y = 0, z = 0 },
                                         new Vector3d { x = 0, y = 1, z = 0 },
                                         new Vector3d { x = 0, y = 0, z = 1 } };
    private NBodySystem system;
    private Dictionary<string, LineRenderer> trajectories;
    private Vector3d Δq;
    private Vector3d Δv;
    private enum ReferenceFrameType {
      [Description("Body-Centric Inertial")]
      Inertial,
      [Description("Body-Centric Rotating")]
      Surface,
      [Description("Barycentric Co-Rotating")]
      CoRotating
    }

    #region Unity Methods

    private void FixedUpdate() {
      if (!FlightGlobals.ready) {
        return;
      }
      Vessel activeVessel = FlightGlobals.ActiveVessel;
      double UT = Planetarium.GetUniversalTime();
      QuaternionD rotation = Planetarium.Rotation;
      if (activeVessel.loaded && !activeVessel.packed) {
        // We compute the velocity ourselves so that we only deal with Unity's
        // and not with KSP's calculations. This also ensures that the
        // accumulation is done in double precision.
        Vector3d newVelocity = Vector3d.zero;
        double totalMass = 0;
        foreach (Part part in activeVessel.parts) {
          newVelocity += (Vector3d)part.rb.velocity * (double)part.rb.mass;
          totalMass += (double)part.rb.mass;
        }
        newVelocity /= totalMass;
        newVelocity += Krakensbane.GetFrameVelocity();

        if (lastFrameWasPhysical) {
          // We now extract the proper acceleration computed by Unity.
          activeVesselProperAcceleration = (newVelocity - activeVesselVelocity)
                              / TimeWarp.fixedDeltaTime - geometricAcceleration;
        } else {
          lastFrameWasPhysical = true;
          activeVesselProperAcceleration = Vector3d.zero;
        }
        // We compute the geometric acceleration which will be applied by KSP.
        CelestialBody primary = activeVessel.orbit.referenceBody;
        Vector3d vesselPosition = activeVessel.findLocalCenterOfMass();
        Vector3d vesselVelocity =
          ((Vector3d)activeVessel.rootPart.rb.GetPointVelocity(vesselPosition))
          + Krakensbane.GetFrameVelocity();
        geometricAcceleration =
          FlightGlobals.getGeeForceAtPosition(vesselPosition)
          + FlightGlobals.getCoriolisAcc(vesselVelocity, primary)
          + FlightGlobals.getCentrifugalAcc(vesselPosition, primary);
        activeVesselVelocity = newVelocity;
      } else {
        lastFrameWasPhysical = false;
        activeVesselProperAcceleration = Vector3d.zero;
      }

      properAccelerationVisualisation.SetPosition(0, activeVessel.rootPart.rb.worldCenterOfMass);
      properAccelerationVisualisation.SetPosition(1, activeVessel.rootPart.rb.worldCenterOfMass + activeVesselProperAcceleration * 100);

      if (simulate) {
#if TRACE
        print("Initiating simulation step...");
#endif
        {
          Body vessel;
          bodies.TryGetValue(activeVessel.id.ToString(), out vessel);
          vessel.properAcceleration = (QuaternionD.Inverse(rotation)
            * activeVesselProperAcceleration).xzy.ToCoordinates();
        }
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
#if TRACE
            print("Updating " + vessel.name + "...");
#endif
          Body secondary, primary, sun;
          bodies.TryGetValue(vessel.id.ToString(), out secondary);
          bodies.TryGetValue(vessel.orbit.referenceBody.name,
                             out primary);
          bodies.TryGetValue("Sun", out sun);
          Vector3d relativePosition
            = secondary.q.ToVector() - primary.q.ToVector();
          Vector3d relativeVelocity
            = secondary.v.ToVector() - primary.v.ToVector();
          if (vessel.isActiveVessel &&
              !vessel.packed) {
            vessel.orbitDriver.UpdateOrbit();
            Vector3d worldRelativePosition = rotation * relativePosition.xzy;
            Vector3d actualRelativePosition = vessel.orbit.getRelativePositionAtUT(UT).xzy;
            Δq = worldRelativePosition - actualRelativePosition;
            vessel.Translate(Δq);
            // We change the velocity ourselves so that we are sure the change
            // is done *now*, and affects the total momentum correctly without
            // changing the relative velocities.
            Vector3d frameVelocity =
              vessel.orbit.referenceBody.inverseRotation ?
                vessel.orbit.referenceBody.getRFrmVel(vessel.orbit.pos) :
                Vector3d.zero;
            Δv = rotation * relativeVelocity.xzy - vessel.orbit.getOrbitalVelocityAtUT(UT).xzy;
            Vector3d Δqplanet =
              rotation * (primary.q.ToVector() - sun.q.ToVector()).xzy
              - vessel.orbit.referenceBody.orbit.getRelativePositionAtUT(UT).xzy;
            print("{" + Δq.ToMathematica() + ","
                      + Δv.ToMathematica() + ","
                      + Δqplanet.ToMathematica() + "}, ");
            foreach (Part part in vessel.parts) {
              part.rb.velocity += Δv;
            }
          } else {
            vessel.orbit.UpdateFromStateVectors(
              (rotation * relativePosition.xzy).xzy,
              (rotation * relativeVelocity.xzy).xzy,
              vessel.orbit.referenceBody, UT);
          }
        }
      }

      if (activeVessel.loaded && !activeVessel.packed) {
        // Compute the velocity used to extract the proper acceleration in the
        // next step.
        activeVesselVelocity = Vector3d.zero;
        double totalMass = 0;
        foreach (Part part in activeVessel.parts) {
          activeVesselVelocity += (Vector3d)part.rb.velocity
            * (double)part.rb.mass;
          totalMass += (double)part.rb.mass;
        }
        activeVesselVelocity /= totalMass;
        activeVesselVelocity += Krakensbane.GetFrameVelocity();
      }

      if (predictTrajectories) {
        system.AdvancePredictions(tmax: UT + 1 * Day,
                                  maxTimestep: 10 * Second,
                                  samplingPeriod: 10);
      }
    }
    private void OnDestroy() {
      RenderingManager.RemoveFromPostDrawQueue(3, new Callback(DrawGUI));
    }
    private void Start() {
      print("Principia: Starting...");
      bodies = new Dictionary<string, Body>();
      trajectories = new Dictionary<string, LineRenderer>();
      referenceBody = Planetarium.fetch.Sun;
      referenceFrameType = ReferenceFrameType.Inertial;
      RenderingManager.AddToPostDrawQueue(3, new Callback(DrawGUI));
      if ((mainWindowPosition.x == 0)
          && (mainWindowPosition.y == 0)) {
        mainWindowPosition = new Rect(Screen.width / 2,
                                      Screen.height / 2,
                                      10, 10);
      }
      if ((referenceFrameWindowPosition.x == 0)
          && (referenceFrameWindowPosition.y == 0)) {
        referenceFrameWindowPosition = new Rect(Screen.width / 3,
                                                Screen.height / 3,
                                                10, 10);
      }
      GameObject gameObject = new GameObject("Line");
      properAccelerationVisualisation = gameObject.AddComponent<LineRenderer>();
      properAccelerationVisualisation.transform.parent = transform;
      properAccelerationVisualisation.useWorldSpace = true;
      properAccelerationVisualisation.transform.localPosition = Vector3.zero;
      properAccelerationVisualisation.transform.localEulerAngles = Vector3.zero;
      properAccelerationVisualisation.material = new Material(Shader.Find("Particles/Additive"));
      properAccelerationVisualisation.SetColors(Color.yellow, Color.green);
      properAccelerationVisualisation.SetWidth(1, 0);
      properAccelerationVisualisation.SetVertexCount(2);

      mapRenderer = MapRenderer.CreateAndAttach(DrawTrajectories);
      print("Principia: Started.");
    }
    private void Update() {
      DrawTrajectories();
    }

    #endregion Unity Methods

    #region GUI

    private void DrawGUI() {
      GUI.skin = HighLogic.Skin;
      mainWindowPosition = GUILayout.Window(1,
                                   mainWindowPosition,
                                   MainWindow,
                                   "Traces of Various Descriptions",
                                   GUILayout.MinWidth(500));
      referenceFrameWindowPosition = GUILayout.Window(2,
                                     referenceFrameWindowPosition,
                                     ReferenceFrameWindow,
                                     "Reference Frame",
                                     GUILayout.MinWidth(500));
    }
    private void MainWindow(int windowID) {
      GUIStyle mySty = new GUIStyle(GUI.skin.button);
      mySty.normal.textColor = mySty.focused.textColor = Color.white;
      mySty.hover.textColor = mySty.active.textColor = Color.yellow;
      mySty.onNormal.textColor
        = mySty.onFocused.textColor
        = mySty.onHover.textColor
        = mySty.onActive.textColor = Color.green;
      mySty.padding = new RectOffset(8, 8, 8, 8);

      GUILayout.BeginVertical();
      GUILayout.TextArea("Planetarium rotation: "
                         + Planetarium.Rotation.w.ToString("F5") + " + "
                         + Planetarium.Rotation.x.ToString("F5") + " i + "
                         + Planetarium.Rotation.y.ToString("F5") + " j + "
                         + Planetarium.Rotation.z.ToString("F5") + " k",
                         GUILayout.ExpandWidth(true));
      if (GUILayout.Button(simulate ? "Switch to Keplerian"
                                    : "Switch to Newtonian",
                           mySty, GUILayout.ExpandWidth(true))) {
        simulate = !simulate;
        double UT = Planetarium.GetUniversalTime();
        bodies.Clear();
        foreach (CelestialBody body in FlightGlobals.Bodies) {
          print("Adding body " + body.name + "...");
          bodies.Add(body.name, body.ToBody());
        }
        foreach (Vessel vessel in FlightGlobals.Vessels) {
          print("Adding vessel " + vessel.name + "...");
          bodies.Add(vessel.id.ToString(), vessel.ToBody());
        }
        system = new NBodySystem(bodies.Values.ToArray(), UT);
      }
      if (GUILayout.Button(
                    predictTrajectories ? "Stop Plotting" : "Plot Trajectories",
                    mySty, GUILayout.ExpandWidth(true))) {
        predictTrajectories = !predictTrajectories;
        double UT = Planetarium.GetUniversalTime();
        if (predictTrajectories) {
          foreach (Vessel vessel in FlightGlobals.Vessels) {
            print("Adding vessel " + vessel.name + " trajectory...");
            GameObject lineObject = new GameObject("Line");
            // TODO(robin): Switch to layer 31, use a mesh, draw the lines in
            // 2D if MapView.Draw3DLines is false (needs a MapViewLine class).
            lineObject.layer = 9;
            LineRenderer line = lineObject.AddComponent<LineRenderer>();
            line.transform.parent = null;
            line.material = MapView.fetch.orbitLinesMaterial;
            line.SetColors(Color.blue, Color.red);
            line.useWorldSpace = true;
            trajectories.Add(vessel.id.ToString(), line);
          }
          system.RecalculatePredictions(tmax: UT + 1 * Day,
                                        maxTimestep: 10 * Second,
                                        samplingPeriod: 10);
        } else {
          foreach (LineRenderer line in trajectories.Values) {
            Destroy(line);
          }
          trajectories.Clear();
        }
      }
      if (simulate) {
        GUILayout.TextArea("Simulating n-body physics...");
        Body kerbin;
        bodies.TryGetValue("Kerbin", out kerbin);
        CelestialBody minmus = FlightGlobals.Bodies[3];
        GUILayout.TextArea("qKerbin: " + kerbin.q.ToVector().ToString("F5"));
        GUILayout.TextArea("vKerbin: " + kerbin.v.ToVector().ToString("F5"));
        GUILayout.TextArea("ωMinmus: " + minmus.angularVelocity.ToString("F5"));
      }
      GUILayout.TextArea("MapView.Draw3DLines: "
                         + MapView.Draw3DLines.ToString());
      GUILayout.TextArea("ActiveVessel.GetWorldPos3D: " +
                     FlightGlobals.ActiveVessel.GetWorldPos3D().ToString("F5"));
      GUILayout.TextArea("Krakensbane.GetFrameVelocity: " +
                         Krakensbane.GetFrameVelocity().ToString("F5"));
      GUILayout.TextArea("FlightGlobals.ActiveVessel.rb_velocity: "
           + ((Vector3d)FlightGlobals.ActiveVessel.rb_velocity).ToString("F5"));
      GUILayout.TextArea("FlightGlobals.ActiveVessel.rootPart.rb.velocity: "
        + ((Vector3d)FlightGlobals.ActiveVessel.rootPart.rb.velocity)
                                                               .ToString("F5"));
      GUILayout.TextArea("FlightGlobals.ActiveVessel.perturbation: "
          + ((Vector3d)FlightGlobals.ActiveVessel.perturbation).ToString("F5"));
      GUILayout.TextArea("FlightGlobals.ActiveVessel.geeForce: "
                          + FlightGlobals.ActiveVessel.geeForce.ToString());
      GUILayout.TextArea("FlightGlobals.ActiveVessel.geeForce_immediate: "
                + FlightGlobals.ActiveVessel.geeForce_immediate.ToString());
      GUILayout.TextArea("Proper Acceleration: "
        + activeVesselProperAcceleration.ToString("F9"));
      GUILayout.TextArea("Δq: " + Δq.ToString("E13"));
      GUILayout.TextArea("Δv: " + Δv.ToString("E13"));
      string integrators = "";
      foreach (Part part in FlightGlobals.ActiveVessel.parts) {
        integrators += part.partName + ": "
          + (part.flightIntegrator == null ? "null; " : "not null; ");
      }
      GUILayout.TextArea(integrators);
      GUILayout.EndVertical();

      GUI.DragWindow(new Rect(left: 0f, top: 0f, width: 10000f, height: 20f));
    }
    private void ReferenceFrameWindow(int windowID) {
      GUIStyle mySty = new GUIStyle(GUI.skin.button);
      mySty.normal.textColor = mySty.focused.textColor = Color.white;
      mySty.hover.textColor = mySty.active.textColor = Color.yellow;
      mySty.onNormal.textColor
        = mySty.onFocused.textColor
        = mySty.onHover.textColor
        = mySty.onActive.textColor = Color.green;
      mySty.padding = new RectOffset(8, 8, 8, 8);

      GUILayout.BeginHorizontal();
      GUILayout.BeginVertical(GUILayout.Width(250));
      foreach (CelestialBody body in FlightGlobals.Bodies) {
        if (GUILayout.Toggle(referenceBody.name == body.name,
                             body.name)
            && referenceBody.name != body.name) {
          referenceBody = body;
        }
      }
      GUILayout.EndVertical();
      GUILayout.BeginVertical(GUILayout.ExpandWidth(true));
      foreach (ReferenceFrameType type in
               Enum.GetValues(typeof(ReferenceFrameType))) {
        if (GUILayout.Toggle(referenceFrameType == type,
                             type.Description(),
                             GUILayout.ExpandWidth(true))
            && referenceFrameType != type
            && referenceBody.name != "Sun") {
          referenceFrameType = type;
        }
      }
      GUILayout.EndVertical();
      GUILayout.EndHorizontal();

      GUI.DragWindow(new Rect(left: 0f, top: 0f, width: 10000f, height: 20f));
    }

    #endregion GUI

    private void DrawTrajectories() {
      if (MapView.MapIsEnabled && predictTrajectories) {
        QuaternionD rotation = Planetarium.Rotation;
        PlanetariumCamera camera = (PlanetariumCamera)
          GameObject.FindObjectOfType(typeof(PlanetariumCamera));
        double UT = Planetarium.GetUniversalTime();
        foreach (Vessel v in FlightGlobals.Vessels) {
          // TODO(robin): There are too many nested expressions here, and I've
          // caught quite a few bugs from misplaced parentheses. Declare more
          // variables so this part is at least vaguely readable.
          LineRenderer line;
          Body vessel;
          Body reference;
          Body primary = null;              // The silly language forbids the
          CelestialBody primaryBody = null; // use of 'unassigned' variables.
          trajectories.TryGetValue(v.id.ToString(), out line);
          bodies.TryGetValue(v.id.ToString(), out vessel);
          bodies.TryGetValue(referenceBody.name, out reference);
          if (referenceFrameType == ReferenceFrameType.CoRotating) {
            primaryBody = referenceBody.orbit.referenceBody;
            bodies.TryGetValue(primaryBody.name, out primary);
          }
          line.SetVertexCount(vessel.predictedTrajectory.Count);
          line.enabled = true;
          Vector3d currentVesselPosition = v.GetWorldPos3D();
          line.SetWidth((0.01f * camera.Distance), (0.01f * camera.Distance));
#pragma warning disable 618 // Obsolete.
          // Quaternion.AxisAngle uses radians (which we want, like and use)
          // and is deprecated because of that for some mindbogglingly stupid
          // reason. We use it nonetheless.
          // We also use Quaternion rather than QuaternionD because KSP's
          // QuaternionD seems broken (maybe it doesn't implement deprecated
          // functions).
          // TODO(robin): Make a versor; I wouldn't want to trust Unity with a
          // calculation even if it were done in double precision.
          for (int i = 0; i < vessel.predictedTrajectory.Count; ++i) {
            switch (referenceFrameType) {
              case ReferenceFrameType.CoRotating:
                Vector3d predictedBarycenter =
                  (reference.predictedTrajectory[i].q.ToVector().xzy
                  * referenceBody.Mass
                  + primary.predictedTrajectory[i].q.ToVector().xzy
                  * primaryBody.Mass)
                / (referenceBody.Mass + primaryBody.Mass);
                Vector3d currentBarycenter =
                  (reference.q.ToVector().xzy * referenceBody.Mass
                  + primary.q.ToVector().xzy * primaryBody.Mass)
                  / (referenceBody.Mass + primaryBody.Mass);
                line.SetPosition(i, ScaledSpace.LocalToScaledSpace(
                  rotation * (
                  Quaternion.FromToRotation(
                    reference.predictedTrajectory[i].q.ToVector().xzy
                      - predictedBarycenter,
                    reference.q.ToVector().xzy - currentBarycenter)
                  * (vessel.predictedTrajectory[i].q.ToVector().xzy
                    - predictedBarycenter) + currentBarycenter
                  - vessel.q.ToVector().xzy) + currentVesselPosition));
                break;
              case ReferenceFrameType.Inertial:
                line.SetPosition(i, ScaledSpace.LocalToScaledSpace(
                  rotation * (vessel.predictedTrajectory[i].q.ToVector().xzy
                           - reference.predictedTrajectory[i].q.ToVector().xzy
                           + reference.q.ToVector().xzy
                           - vessel.q.ToVector().xzy)
                  + currentVesselPosition));
                break;
              case ReferenceFrameType.Surface:
                line.SetPosition(i, ScaledSpace.LocalToScaledSpace(
                  rotation * (
                  Quaternion.AxisAngle(referenceBody.angularVelocity, (float)(
                                       referenceBody.angularV * (UT -
                                            vessel.predictedTrajectory[i].t)))
                    * (vessel.predictedTrajectory[i].q.ToVector().xzy
                      - reference.predictedTrajectory[i].q.ToVector().xzy)
                    + reference.q.ToVector().xzy - vessel.q.ToVector().xzy)
                  + currentVesselPosition));
                break;
            }
          }
#pragma warning restore 618
        }
      }
    }
  }
}