/*
 * Copyright� (c) 2018 Maarten Maathuis, (aka madman2003).
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.

 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

using System;
using System.Collections.Generic;
using System.Linq;

namespace principia {
namespace ksp_plugin_adapter {

    static class CelestialExtensions {
        public static bool is_leaf(this CelestialBody celestial) {
            return celestial.orbitingBodies.Count == 0;
        }

        public static bool is_root(this CelestialBody celestial) {
            return celestial.orbit == null;
        }
    }

    // For the sake of not having to know the GUI complexity when working with the controllers
    // that actually do stuff with data, we isolate the data into a seperate class
    // It's intended to be as simple as possible, so that at any time both the GUI and control
    // developers can easily read this
    public sealed class DataServices
    {
        //
        // Global settings
        //
        private static IntPtr plugin;
        private static ManeuverNode guidance_node;

        public static void SetPlugin(IntPtr value)
        {
            plugin = value;
            InitializeLoggingSettings();
        }

        private static Vessel GetVessel()
        {
            return FlightGlobals.ActiveVessel;
        }

        private static string GetVesselGuid()
        {
            Vessel vessel = GetVessel();
            if (vessel)
                return GetVessel().id.ToString();
            return "";
        }

        private static bool HasFlightPlan()
        {
            if (!GetVessel())
                return false;

            string vesselguid = GetVesselGuid();
            return plugin.HasVessel(vesselguid) && plugin.FlightPlanExists(vesselguid);
        }

        //
        // History length
        //
        private static int history_magnitude = 20; // (1 << index) is the history time in seconds, with the exception of 30, which is +infinity
        public static float GetHistoryMagnitude() { return (float)history_magnitude; }
        public static void SetHistoryMagnitude(float value) { history_magnitude = (int)value; }
        public static double GetHistoryLength()
        {
            if (history_magnitude == 30) {
                return double.PositiveInfinity;
            } else {
                return (1 << history_magnitude);
            }
        }

        //
        // Prediction
        //
        private static int prediction_tolerance_magnitude = -2;
        private static int prediction_step_magnitude = 8;
        public static float GetPredictionToleranceMagnitude() { return (float)prediction_tolerance_magnitude; }
        public static void SetPredictionToleranceMagnitude(float value) { prediction_tolerance_magnitude = (int)value; }
        public static double GetPredictionTolerance() { return Math.Pow(10, prediction_tolerance_magnitude); }
        public static float GetPredictionStepMagnitude() { return (float)prediction_step_magnitude; }
        public static void SetPredictionStepMagnitude(float value) { prediction_step_magnitude = (int)value; }
        public static int GetPredictionStep() { return (1 << prediction_step_magnitude); }

        //
        // KSP settings
        //
        private static bool display_patched_conics = false;
        public static bool GetPatchedConicsEnabled() { return display_patched_conics; }
        public static void SetPatchedConicsEnabled(bool value) { display_patched_conics = value; }
        // For legacy reasons this setting is directly configured on an KSP object
        public static bool GetSolarFlareEnabled() { return Sun.Instance.sunFlare.enabled; }
        public static void SetSolarFlareEnabled(bool value) { Sun.Instance.sunFlare.enabled = value; }

        //
        // Reference frame settings
        //
        private static CelestialBody selected_celestial_body = FlightGlobals.GetHomeBody();
        // I suspect the hardcoded numbers are to allow the C++ code to understand this as well
        public enum FrameType {
            BARYCENTRIC_ROTATING = 6001,
            BODY_CENTRED_NON_ROTATING = 6000,
            BODY_CENTRED_PARENT_DIRECTION = 6002,
            BODY_SURFACE = 6003
        }
        private static FrameType reference_frame = FrameType.BODY_CENTRED_NON_ROTATING;
        private static FrameType last_non_surface_reference_frame = reference_frame;
        private static Vessel reference_frame_target_override = null;
        // internal is inherited from NavigationFrameParameters
        internal delegate void NavigationFrameParametersCallback(NavigationFrameParameters frame_parameters);
        private static NavigationFrameParametersCallback on_update_celestial_body_or_reference_frame;

        public static CelestialBody GetSelectedCelestialBody() { return selected_celestial_body; }
        public static void SetSelectedCelestialBody(CelestialBody value)
        {
            CelestialBody prev = selected_celestial_body;
            selected_celestial_body = value;
            if (prev != selected_celestial_body)
            {
                UpdatePluginWithCelestialBodyAndReferenceFrame();
            }
        }
        public static FrameType GetReferenceFrame() { return reference_frame; }
        public static void SetReferenceFrame(FrameType value)
        {
            FrameType prev = reference_frame;
            reference_frame = value;
            // We cannot support frames that require 2 bodies, if we selected the root body
            // (and thus only have 1 body to work with)
            if (selected_celestial_body.is_root() &&
               (reference_frame == FrameType.BARYCENTRIC_ROTATING ||
                reference_frame == FrameType.BODY_CENTRED_PARENT_DIRECTION)) {
                reference_frame = FrameType.BODY_CENTRED_NON_ROTATING;
            }
            if (reference_frame != FrameType.BODY_SURFACE)
            {
                last_non_surface_reference_frame = reference_frame;
            }
            if (prev != reference_frame)
            {
                UpdatePluginWithCelestialBodyAndReferenceFrame();
            }
        }
        public static void SetLastNonSurfaceReferenceFrame()
        {
            SetReferenceFrame(last_non_surface_reference_frame);
        }
        public static Vessel GetReferenceFrameTargetOverride()
        {
            return reference_frame_target_override;
        }
        public static void SetReferenceFrameTargetOverride(Vessel value)
        {
            reference_frame_target_override = value;
        }

        // internal is inherited from NavigationFrameParameters
        internal static void InitializeSelectedCelestialBodyAndReferenceFrame(NavigationFrameParametersCallback callback)
        {
            on_update_celestial_body_or_reference_frame = callback;
            reference_frame = FrameType.BODY_CENTRED_NON_ROTATING;
            selected_celestial_body = FlightGlobals.currentMainBody ?? FlightGlobals.GetHomeBody();
        }
        
        private static NavigationFrameParameters GenerateNavigationFrameParameters()
        {
            switch (reference_frame) {
                case FrameType.BODY_CENTRED_NON_ROTATING:
                case FrameType.BODY_SURFACE:
                    return new NavigationFrameParameters{
                        extension = (int)reference_frame,
                        centre_index = selected_celestial_body.flightGlobalsIndex};
                case FrameType.BARYCENTRIC_ROTATING:
                    return new NavigationFrameParameters{
                        extension = (int)reference_frame,
                        primary_index = selected_celestial_body.referenceBody.flightGlobalsIndex,
                        secondary_index = selected_celestial_body.flightGlobalsIndex};
                case FrameType.BODY_CENTRED_PARENT_DIRECTION:
                    // We put the primary body as secondary, because the one we want fixed
                    // is the secondary body (which means it has to be the primary in the
                    // terminology of |BodyCentredBodyDirection|).
                    return new NavigationFrameParameters{
                        extension = (int)reference_frame,
                        primary_index = selected_celestial_body.flightGlobalsIndex,
                        secondary_index = selected_celestial_body.referenceBody.flightGlobalsIndex};
                default:
                    throw Log.Fatal("Unexpected reference_frame " + reference_frame.ToString());
            }
        }

        public static String ReferenceFrameShortName() {
            CelestialBody selected = GetSelectedCelestialBody();
            if (GetReferenceFrameTargetOverride()) {
                return "Tgt LVLH@" + selected.name[0];
            }
            switch (GetReferenceFrame()) {
                case FrameType.BODY_CENTRED_NON_ROTATING:
                    return selected.name[0] + "CI";
                case FrameType.BARYCENTRIC_ROTATING:
                    if (selected.is_root()) {
                        throw Log.Fatal("Naming barycentric rotating frame of root body");
                    } else {
                        return selected.referenceBody.name[0] + (selected.name[0] + "B");
                    }
                case FrameType.BODY_CENTRED_PARENT_DIRECTION:
                    if (selected.is_root()) {
                        throw Log.Fatal("Naming parent-direction rotating frame of root body");
                    } else {
                        return selected.name[0] + "C" + selected.referenceBody.name[0] + "A";
                    }
                case FrameType.BODY_SURFACE:
                    return selected.name[0] + "C" + selected.name[0] + "F";
                default:
                    throw Log.Fatal("Unexpected type " + GetReferenceFrame().ToString());
            }
        }

        public static CelestialBody[] ReferenceFrameFixedBodies() {
            CelestialBody selected = GetSelectedCelestialBody();
            if (GetReferenceFrameTargetOverride()) {
                return new CelestialBody[]{};
            }
            switch (GetReferenceFrame()) {
                case FrameType.BODY_CENTRED_NON_ROTATING:
                case FrameType.BODY_CENTRED_PARENT_DIRECTION:
                case FrameType.BODY_SURFACE:
                    return new CelestialBody[]{selected};
                case FrameType.BARYCENTRIC_ROTATING:
                    return new CelestialBody[]{};
                default:
                    throw Log.Fatal("Unexpected frame_type " + GetReferenceFrame().ToString());
            }
        }

        private static void UpdatePluginWithCelestialBodyAndReferenceFrame()
        {
            NavigationFrameParameters frame_parameters = GenerateNavigationFrameParameters();
            on_update_celestial_body_or_reference_frame(frame_parameters);
        }

        //
        // Logging settings
        //
        private static bool record_journal_at_next_startup = false;
        private static bool record_journal_in_progress = false;
        private static int verbose_level = 0;
        private static int supressed_logging_level = 0;
        private static int stderr_logging_level = 2;
        private static int flush_logging_level = -1;

        public static bool GetRecordJournalInProgress() { return record_journal_in_progress; }
        public static bool GetRecordJournalAtNextStartup() { return record_journal_at_next_startup; }
        public static void SetRecordJournalAtNextStartup(bool value) { record_journal_at_next_startup = value; }
        
        public static void InitializeJournaling()
        {
            if (record_journal_at_next_startup) {
                Log.ActivateRecorder(true);
                record_journal_in_progress = true;
            }
        }

        // Implicit requirement: Don't accept a log setting, until the C++ side of principia has accepted it
        public static int GetVerboseLevel() { return verbose_level; }
        public static void SetVerboseLevel(int value) { Log.SetVerboseLogging(value); verbose_level = Log.GetVerboseLogging(); }
        public static int GetLogLevel() { return supressed_logging_level; }
        public static void SetLogLevel(int value) { Log.SetSuppressedLogging(value); supressed_logging_level = Log.GetSuppressedLogging(); }
        public static int GetStderrLevel() { return stderr_logging_level; }
        public static void SetStderrLevel(int value) { Log.SetStderrLogging(value); stderr_logging_level = Log.GetStderrLogging(); }
        public static int GetFlushLevel() { return flush_logging_level; }
        public static void SetFlushLevel(int value) { Log.SetBufferedLogging(value); flush_logging_level = Log.GetBufferedLogging(); }

        private static void InitializeLoggingSettings()
        {
            SetVerboseLevel(verbose_level);
            SetLogLevel(supressed_logging_level);
            SetStderrLevel(stderr_logging_level);
            SetFlushLevel(flush_logging_level);
        }

        //
        // Planning
        //

        //
        // Planner: settings
        //
        public const double PLAN_TIME_LENGTH_DEFAULT = 10*60; // TODO: is this the proper default?
        public const int PLAN_MAX_STEPS_PER_SEGMENT_DEFAULT = 1000;
        public const double PLAN_TOLERANCE_DEFAULT = 1.0;

        private static double plan_time_length = PLAN_TIME_LENGTH_DEFAULT;
        private static int plan_max_steps_per_segment = PLAN_MAX_STEPS_PER_SEGMENT_DEFAULT;
        private static double plan_tolerance = PLAN_TOLERANCE_DEFAULT;

        public static double GetPlanTimeLength() { return plan_time_length; }
        public static void SetPlanTimeLength(double value)
        {
            plan_time_length = value;
            if (plan_time_length < 10.0) plan_time_length = 10.0;
            UpdateFlightPlanTimeLength();
        }
        public static int GetPlanMaxStepsPerSegment() { return plan_max_steps_per_segment; }
        public static void SetPlanMaxStepsPerSegment(int value)
        {
            plan_max_steps_per_segment = value;
            if (plan_max_steps_per_segment < 100) plan_max_steps_per_segment = 100;
            UpdateFlightPlanAdaptiveStepParameters();
        }
        public static double GetPlanTolerance() { return plan_tolerance; }
        public static void SetPlanTolerance(double value)
        {
            plan_tolerance = value;
            if (plan_tolerance < 1e-6) plan_tolerance = 1e-6;
            if (plan_tolerance > 1e6) plan_tolerance = 1e6;
            UpdateFlightPlanAdaptiveStepParameters();
        }

        private static void UpdateFlightPlanAdaptiveStepParameters()
        {
            if (HasFlightPlan())
            {
                string vesselguid = GetVesselGuid();
                FlightPlanAdaptiveStepParameters parameters = plugin.FlightPlanGetAdaptiveStepParameters(vesselguid);
                parameters.length_integration_tolerance = plan_tolerance;
                parameters.speed_integration_tolerance = plan_tolerance;
                parameters.max_steps = plan_max_steps_per_segment;
                plugin.FlightPlanSetAdaptiveStepParameters(vesselguid, parameters);
            }
        }

        private static void UpdateFlightPlanTimeLength()
        {
            if (HasFlightPlan())
            {
                string vesselguid = GetVesselGuid();
                // TODO: there is also an actual final time, what to do with that beast? some feedback in the GUI needed?
                // TODO: how frequently do we need to set this again?
                plugin.FlightPlanSetDesiredFinalTime(vesselguid, plugin.CurrentTime() + plan_time_length);
                plan_time_length = plugin.FlightPlanGetDesiredFinalTime(vesselguid) - plugin.CurrentTime();
            }
        }

        public static void EnsureFlightPlanExists()
        {
            if (!HasFlightPlan())
            {
                plugin.FlightPlanCreate(GetVesselGuid(),
                                        plugin.CurrentTime() + plan_time_length,
                                        GetVessel().GetTotalMass());
                UpdateFlightPlanAdaptiveStepParameters();
                UpdateFlightPlanTimeLength();
            }
        }

        public static void UpdateSettingsUponVesselChange()
        {
            if (HasFlightPlan())
            {
                UpdateFlightPlanTimeLength();
                UpdateFlightPlanAdaptiveStepParameters();
            }
        }

        //
        // Planner: execution
        //
        private static bool show_on_navball = false;

        public static bool GetShowOnNavball() { return show_on_navball; }
        public static void SetShowOnNavball(bool value) { show_on_navball = value; if (show_on_navball) ShowOnNavball(); else HideFromNavball(); }
        public static double GetEngineDeltaTime()
        {
            int selected_maneuver;
            NavigationManoeuvre maneuver;
            if (!FindUpcomingManeuver(out selected_maneuver, out maneuver))
                return 0.0;

            if (plugin.CurrentTime() < maneuver.burn.initial_time)
                return plugin.CurrentTime() - maneuver.burn.initial_time;
            else
                return plugin.CurrentTime() - maneuver.final_time;
        }
        public static bool IsEngineDeltaTimeIgnition()
        {
            int selected_maneuver;
            NavigationManoeuvre maneuver;
            // People who see a fixed time of 0, are not likely to ignite engines
            // and pretending by default to show engine cutoff time can also be confusing
            if (!FindUpcomingManeuver(out selected_maneuver, out maneuver))
                return true;

            if (plugin.CurrentTime() < maneuver.burn.initial_time)
                return true;
            else
                return false;
        }

        public static double GetDeltaVelocityOfAllBurns()
        {
            double total_delta_v = 0.0;

            if (HasFlightPlan())
            {
                string vesselguid = GetVesselGuid();
                for (int index = 0; index < plugin.FlightPlanNumberOfManoeuvres(vesselguid); index++)
                {
                    total_delta_v += ((Vector3d)plugin.FlightPlanGetManoeuvre(vesselguid, index).burn.delta_v).magnitude;
                }
            }

            return total_delta_v;
        }

        private static bool FindUpcomingManeuver(out int index, out NavigationManoeuvre maneuver)
        {
            index = 0;
            maneuver = new NavigationManoeuvre();

            if (HasFlightPlan())
            {
                string vesselguid = GetVesselGuid();
                int number_of_maneuvers = plugin.FlightPlanNumberOfManoeuvres(vesselguid);

                if (number_of_maneuvers < 1)
                    return false;

                maneuver = plugin.FlightPlanGetManoeuvre(vesselguid, index);
                while ((plugin.CurrentTime() > maneuver.final_time) && ((index + 1) < number_of_maneuvers))
                {
                    index += 1;
                    maneuver = plugin.FlightPlanGetManoeuvre(vesselguid, index);
                }

                // No maneuver satifies the condition
                if (index == number_of_maneuvers)
                    return false;
                return true;
            }
            return false;
        }

        private static bool PossibleToShowOnNavball(Vessel vessel)
        {
            return vessel && (vessel.patchedConicSolver != null);
        }

        private static void ShowOnNavball()
        {
            Vessel vessel = GetVessel();
            int selected_maneuver;
            NavigationManoeuvre maneuver;

            if (!FindUpcomingManeuver(out selected_maneuver, out maneuver))
            {
                HideFromNavball();
                return;
            }
            string vesselguid = GetVesselGuid();

            // In career mode, the patched conic solver may be null.  In that case
            // we do not offer the option of showing the manoeuvre on the navball,
            // even though the flight planner is still available to plan it.
            // TODO(egg): We may want to consider setting the burn vector directly
            // rather than going through the solver.
            if (PossibleToShowOnNavball(vessel)) {
                XYZ guidance = plugin.FlightPlanGetGuidance(vesselguid, selected_maneuver);
                if (!double.IsNaN(guidance.x + guidance.y + guidance.z)) {
                    if (guidance_node == null ||
                        !vessel.patchedConicSolver.maneuverNodes.Contains(guidance_node)) {
                        while (vessel.patchedConicSolver.maneuverNodes.Count > 0) {
                            vessel.patchedConicSolver.maneuverNodes.Last().RemoveSelf();
                        }
                        guidance_node = vessel.patchedConicSolver.AddManeuverNode(maneuver.burn.initial_time);
                    } else if (vessel.patchedConicSolver.maneuverNodes.Count > 1) {
                        while (vessel.patchedConicSolver.maneuverNodes.Count > 1) {
                            if (vessel.patchedConicSolver.maneuverNodes.First() == guidance_node) {
                                vessel.patchedConicSolver.maneuverNodes.Last().RemoveSelf();
                            } else {
                                vessel.patchedConicSolver.maneuverNodes.First().RemoveSelf();
                            }
                        }
                    }
                }

                var stock_orbit = guidance_node.patch;
                Vector3d stock_velocity_at_node_time =
                    stock_orbit.getOrbitalVelocityAtUT(maneuver.burn.initial_time).xzy;
                Vector3d stock_displacement_from_parent_at_node_time =
                    stock_orbit.getRelativePositionAtUT(maneuver.burn.initial_time).xzy;
                UnityEngine.Quaternion stock_frenet_frame_to_world =
                    UnityEngine.Quaternion.LookRotation(
                        stock_velocity_at_node_time,
                        Vector3d.Cross(stock_velocity_at_node_time, stock_displacement_from_parent_at_node_time)
                    );
                guidance_node.DeltaV =
                    ((Vector3d)maneuver.burn.delta_v).magnitude *
                     (Vector3d)(UnityEngine.Quaternion.Inverse(stock_frenet_frame_to_world) *
                     (Vector3d)guidance);
                guidance_node.UT = maneuver.burn.initial_time;
                vessel.patchedConicSolver.UpdateFlightPlan();
            }
        }

        private static void HideFromNavball()
        {
            if (guidance_node != null) {
                guidance_node.RemoveSelf();
                guidance_node = null;
            }
        }

        private static void RefreshShowOnNavball()
        {
            if (show_on_navball)
            {
                ShowOnNavball(); // refresh
            }
        }

        public static void WarpToManeuver()
        {
            int selected_maneuver;
            NavigationManoeuvre maneuver;
            if (!FindUpcomingManeuver(out selected_maneuver, out maneuver))
                return;
            TimeWarp.fetch.WarpTo(maneuver.burn.initial_time - 60);
        }

        //
        // Planner: plan
        //
        public static double DEFAULT_MANEUVER_DELTA_TIME = 60.0;
        public enum BurnMode {Engine, RCS, Instant};
        private static List<BurnMode> burn_mode = new List<BurnMode>();

        private static NavigationManoeuvre GetFlightPlanManoeuver(int index)
        {
            if (HasFlightPlan() && GetLastManeuverIndex() >= index)
            {
                return plugin.FlightPlanGetManoeuvre(GetVesselGuid(), index);
            }
            return new NavigationManoeuvre(); // empty placeholder until we can provide a meaningful non-zero value
        }

        public static double GetManeuverDeltaVelocityTangent(int index) { return GetFlightPlanManoeuver(index).burn.delta_v.x; }
        public static void SetManeuverDeltaVelocityTangent(double value)
        {
            int last_index = GetLastManeuverIndex();
            if (last_index < 0)
                return;
            NavigationManoeuvre last = GetFlightPlanManoeuver(last_index);
            last.burn.delta_v.x = value;
            plugin.FlightPlanReplaceLast(GetVesselGuid(), last.burn);
            RefreshShowOnNavball();
        }
        public static double GetManeuverDeltaVelocityNormal(int index) { return GetFlightPlanManoeuver(index).burn.delta_v.y; }
        public static void SetManeuverDeltaVelocityNormal(double value)
        {
            int last_index = GetLastManeuverIndex();
            if (last_index < 0)
                return;
            NavigationManoeuvre last = GetFlightPlanManoeuver(last_index);
            last.burn.delta_v.y = value;
            plugin.FlightPlanReplaceLast(GetVesselGuid(), last.burn);
            RefreshShowOnNavball();
        }
        public static double GetManeuverDeltaVelocityBinormal(int index) { return GetFlightPlanManoeuver(index).burn.delta_v.z; }
        public static void SetManeuverDeltaVelocityBinormal(double value)
        {
            int last_index = GetLastManeuverIndex();
            if (last_index < 0)
                return;
            NavigationManoeuvre last = GetFlightPlanManoeuver(last_index);
            last.burn.delta_v.z = value;
            plugin.FlightPlanReplaceLast(GetVesselGuid(), last.burn);
            RefreshShowOnNavball();
        }
        public static double GetManeuverDeltaTime(int index) { return GetFlightPlanManoeuver(index).burn.initial_time - plugin.CurrentTime(); }
        public static void SetManeuverDeltaTime(double value)
        {
            int last_index = GetLastManeuverIndex();
            if (last_index < 0)
                return;
            NavigationManoeuvre last = GetFlightPlanManoeuver(last_index);
            last.burn.initial_time = value + plugin.CurrentTime();
            plugin.FlightPlanReplaceLast(GetVesselGuid(), last.burn);
            RefreshShowOnNavball();
        }
        public static double GetManeuverTime(int index) { return GetFlightPlanManoeuver(index).burn.initial_time; }
        public static bool GetManeuverIntertiallyFixed(int index) { return GetFlightPlanManoeuver(index).burn.is_inertially_fixed; }
        public static void SetManeuverIntertiallyFixed(bool value)
        {
            int last_index = GetLastManeuverIndex();
            if (last_index < 0)
                return;
            NavigationManoeuvre last = GetFlightPlanManoeuver(last_index);
            last.burn.is_inertially_fixed = value;
            plugin.FlightPlanReplaceLast(GetVesselGuid(), last.burn);
            RefreshShowOnNavball();
        }

        public static double GetBurnDeltaVelocity(int index)
        {
            NavigationManoeuvre maneuver = GetFlightPlanManoeuver(index);
            return new Vector3d{x = maneuver.burn.delta_v.x,
                                y = maneuver.burn.delta_v.y,
                                z = maneuver.burn.delta_v.z}.magnitude;
        }
        public static double GetBurnTime(int index) { return GetFlightPlanManoeuver(index).duration; }

        public static BurnMode GetBurnMode(int index)
        {
            if (index < burn_mode.Count)
                return burn_mode[index];
            return BurnMode.Engine;
        }
        public static void SetBurnMode(BurnMode value)
        {
            int index = GetLastManeuverIndex();
            while (index >= burn_mode.Count) {
                burn_mode.Add(BurnMode.Engine);
            }
            burn_mode[index] = value;

            if (index < 0)
                return;
            NavigationManoeuvre last = GetFlightPlanManoeuver(index);
            EnginePerformance engine_performance = GetEnginePerformance(value);
            last.burn.thrust_in_kilonewtons = engine_performance.thrust_in_kilonewtons;
            last.burn.specific_impulse_in_seconds_g0 = engine_performance.specific_impulse_in_seconds_g0;
            plugin.FlightPlanReplaceLast(GetVesselGuid(), last.burn);

            RefreshShowOnNavball();
        }

        public static int GetLastManeuverIndex()
        {
            if (!HasFlightPlan())
                return -1;
            string vesselguid = GetVesselGuid();
            return plugin.FlightPlanNumberOfManoeuvres(vesselguid) - 1;
        }

        private static EnginePerformance GetEnginePerformance(BurnMode burn_mode)
        {
            EnginePerformance engine_performance;
            switch (burn_mode)
            {
                case BurnMode.Engine:
                    engine_performance = ComputeEngineCharacteristics();
                    break;
                case BurnMode.RCS:
                    engine_performance = ComputeRCSCharacteristics();
                    break;
                case BurnMode.Instant:
                default:
                    engine_performance = UseTheForceLuke();
                    break;
            }
            return engine_performance;
        }

        public static bool AddManeuver()
        {
            if (GetVessel())
            {
                string vesselguid = GetVesselGuid();
                int last_maneuver_index = GetLastManeuverIndex();
                double initial_time = 0.0;
                bool inertially_fixed = false;
                EnginePerformance engine_performance;
                BurnMode burn_mode = BurnMode.Engine;
                bool return_value = false;

                EnsureFlightPlanExists();

                if (last_maneuver_index >= 0)
                {
                    initial_time = Math.Max(plugin.CurrentTime(), GetManeuverTime(last_maneuver_index)) + 60;
                    inertially_fixed = GetManeuverIntertiallyFixed(last_maneuver_index);
                    burn_mode = GetBurnMode(last_maneuver_index);
                }
                else
                {
                    initial_time = plugin.CurrentTime() + 60;
                }

                engine_performance = GetEnginePerformance(burn_mode);
                Burn candidate_burn = new Burn{
                    thrust_in_kilonewtons = engine_performance.thrust_in_kilonewtons,
                    specific_impulse_in_seconds_g0 = engine_performance.specific_impulse_in_seconds_g0,
                    frame = GenerateNavigationFrameParameters(),
                    initial_time = initial_time,
                    delta_v = new XYZ{x = 0.0,
                                      y = 0.0,
                                      z = 0.0},
                    is_inertially_fixed = inertially_fixed};
                return_value = plugin.FlightPlanAppend(vesselguid, candidate_burn);
                // This setting is not derived from the maneuver node, so we must store it
                SetBurnMode(burn_mode);
                RefreshShowOnNavball();
                return return_value;
            }
            return false;
        }

        public static void RemoveLastManeuver()
        {
            if (GetVessel())
            {
                string vesselguid = GetVesselGuid();
                if (GetLastManeuverIndex() >= 0)
                {
                    plugin.FlightPlanRemoveLast(vesselguid);
                }
                if (GetLastManeuverIndex() < 0)
                {
                    plugin.FlightPlanDelete(vesselguid);
                }
                RefreshShowOnNavball();
            }
        }

        //
        // Engine performance calculations
        //

        private class EnginePerformance
        {
            public double thrust_in_kilonewtons;
            public double specific_impulse_in_seconds_g0;
        }

        // TODO: if how should we inform users of fallback scenarios?
        private static EnginePerformance ComputeEngineCharacteristics() {
            EnginePerformance engine_performance = new EnginePerformance();
            Vessel vessel = GetVessel();

            ModuleEngines[] active_engines =
                (from part in vessel.parts
                select (from PartModule module in part.Modules
                        where module is ModuleEngines &&
                              (module as ModuleEngines).EngineIgnited
                        select module as ModuleEngines)).SelectMany(x => x).ToArray();

            Vector3d reference_direction = vessel.ReferenceTransform.up;
            double[] thrusts =
                (from engine in active_engines
                 select engine.maxThrust *
                    (from transform in engine.thrustTransforms
                     select Math.Max(0,
                                     Vector3d.Dot(reference_direction,
                                                  -transform.forward))).Average()).ToArray();
            engine_performance.thrust_in_kilonewtons = thrusts.Sum();

            // This would use zip if we had 4.0 or later.  We loop for now.
            double Σ_f_over_i_sp = 0;
            for (int i = 0; i < active_engines.Count(); ++i) {
                Σ_f_over_i_sp += thrusts[i] / active_engines[i].atmosphereCurve.Evaluate(0);
            }
            engine_performance.specific_impulse_in_seconds_g0 = engine_performance.thrust_in_kilonewtons / Σ_f_over_i_sp;

            // If there are no engines, fall back onto RCS.
            if (engine_performance.thrust_in_kilonewtons == 0) {
                return ComputeRCSCharacteristics();
            }

            return engine_performance;
        }

        private static EnginePerformance ComputeRCSCharacteristics() {
            EnginePerformance engine_performance = new EnginePerformance();
            Vessel vessel = GetVessel();

            ModuleRCS[] active_rcs =
                (from part in vessel.parts
                 select (from PartModule module in part.Modules
                        where module is ModuleRCS &&
                              (module as ModuleRCS).rcsEnabled
                        select module as ModuleRCS)).SelectMany(x => x).ToArray();

            Vector3d reference_direction = vessel.ReferenceTransform.up;
            // NOTE(egg): NathanKell informs me that in >= 1.0.5, RCS has a useZaxis
            // property, that controls whether they thrust -up or -forward.  The madness
            // keeps piling up.
            double[] thrusts =
                (from engine in active_rcs
                 select engine.thrusterPower *
                    (from transform in engine.thrusterTransforms
                     select Math.Max(0,
                                     Vector3d.Dot(reference_direction,
                                                  -transform.up))).Average()).ToArray();
            engine_performance.thrust_in_kilonewtons = thrusts.Sum();

            // This would use zip if we had 4.0 or later.  We loop for now.
            double Σ_f_over_i_sp = 0;
            for (int i = 0; i < active_rcs.Count(); ++i) {
                Σ_f_over_i_sp += thrusts[i] / active_rcs[i].atmosphereCurve.Evaluate(0);
            }
            engine_performance.specific_impulse_in_seconds_g0 = engine_performance.thrust_in_kilonewtons / Σ_f_over_i_sp;

            // If RCS provides no thrust, model a virtually instant burn.
            if (engine_performance.thrust_in_kilonewtons == 0) {
                return UseTheForceLuke();
            }

            return engine_performance;
        }

        private static EnginePerformance UseTheForceLuke() {
            EnginePerformance engine_performance = new EnginePerformance();
            Vessel vessel = GetVessel();

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
            engine_performance.thrust_in_kilonewtons = vessel.GetTotalMass() * range * scale;
            engine_performance.specific_impulse_in_seconds_g0 = range;

            return engine_performance;
        }
    }
}  // namespace ksp_plugin_adapter
}  // namespace principia