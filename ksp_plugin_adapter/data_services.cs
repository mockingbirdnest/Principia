/*
 * Copyright© (c) 2018 Maarten Maathuis, (aka madman2003).
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
using System.Linq;

namespace principia {
namespace ksp_plugin_adapter {

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
        private static bool events_registered = false;
        private static ManeuverNode guidance_node;

        public static void SetPlugin(IntPtr value)
        {
            plugin = value;
            if (!events_registered)
            {
                GameEvents.onVesselChange.Add(OnVesselChange);
                events_registered = true;
            }
            InitializeLoggingSettings();
        }

        private static Vessel GetVessel()
        {
            return FlightGlobals.ActiveVessel;
        }

        private static string GetVesselGuid()
        {
            return GetVessel().id.ToString();
        }

        private static bool HasFlightPlan()
        {
            return plugin.FlightPlanExists(GetVesselGuid());
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
            if (prev != reference_frame)
            {
                UpdatePluginWithCelestialBodyAndReferenceFrame();
            }
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
            string vesselguid = GetVesselGuid();
            if (HasFlightPlan())
            {
                FlightPlanAdaptiveStepParameters parameters = plugin.FlightPlanGetAdaptiveStepParameters(vesselguid);
                parameters.length_integration_tolerance = plan_tolerance;
                parameters.speed_integration_tolerance = plan_tolerance;
                parameters.max_steps = plan_max_steps_per_segment;
                plugin.FlightPlanSetAdaptiveStepParameters(vesselguid, parameters);
            }
        }

        private static void UpdateFlightPlanTimeLength()
        {
            string vesselguid = GetVesselGuid();
            if (HasFlightPlan())
            {
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

        private static void OnVesselChange(Vessel value)
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
            if (!FindUpcomingManeuver(out selected_maneuver, out maneuver))
                return true; // TODO: what should the default be?

            if (plugin.CurrentTime() < maneuver.burn.initial_time)
                return true;
            else
                return false;
        }

        public static double GetDeltaVelocityOfAllBurns()
        {
            string vesselguid = GetVesselGuid();
            double total_delta_v = 0.0;

            if (HasFlightPlan())
            {
                for (int index = 0; index < plugin.FlightPlanNumberOfManoeuvres(vesselguid); index++)
                {
                    total_delta_v += ((Vector3d)plugin.FlightPlanGetManoeuvre(vesselguid, index).burn.delta_v).magnitude;
                }
            }

            return total_delta_v;
        }

        private static bool FindUpcomingManeuver(out int index, out NavigationManoeuvre maneuver)
        {
            string vesselguid = GetVesselGuid();
            index = 0;
            maneuver = new NavigationManoeuvre();

            if (HasFlightPlan())
            {
                int number_of_maneuvers = plugin.FlightPlanNumberOfManoeuvres(vesselguid);

                if (number_of_maneuvers == 0)
                    return false;

                maneuver = plugin.FlightPlanGetManoeuvre(vesselguid, index);
                while ((plugin.CurrentTime() > maneuver.final_time) && (index < number_of_maneuvers))
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

        // TODO: What are the events to update this beast?
        private static void ShowOnNavball()
        {
            Vessel vessel = GetVessel();
            string vesselguid = GetVesselGuid();
            int selected_maneuver;
            NavigationManoeuvre maneuver;

            if (!FindUpcomingManeuver(out selected_maneuver, out maneuver))
                return;

            // In career mode, the patched conic solver may be null.  In that case
            // we do not offer the option of showing the manoeuvre on the navball,
            // even though the flight planner is still available to plan it.
            // TODO(egg): We may want to consider setting the burn vector directly
            // rather than going through the solver.
            if (vessel.patchedConicSolver != null) {
                XYZ guidance = plugin.FlightPlanGetGuidance(vesselguid, selected_maneuver);
                if (!double.IsNaN(guidance.x + guidance.y + guidance.z)) {
                    if (guidance_node == null ||
                        !vessel.patchedConicSolver.maneuverNodes.Contains(guidance_node)) {
                        while (vessel.patchedConicSolver.maneuverNodes.Count > 0) {
                            vessel.patchedConicSolver.maneuverNodes.Last().RemoveSelf();
                        }
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
        public enum BurnMode {Engine, RCS, Instant};
        private static BurnMode burn_mode = BurnMode.Engine;

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
        }

        public static double GetBurnDeltaVelocity(int index)
        {
            NavigationManoeuvre maneuver = GetFlightPlanManoeuver(index);
            return new Vector3d{x = maneuver.burn.delta_v.x,
                                y = maneuver.burn.delta_v.y,
                                z = maneuver.burn.delta_v.z}.magnitude;
        }
        public static double GetBurnTime(int index) { return GetFlightPlanManoeuver(index).duration; }

        // TODO: this needs to be hooked up to engine calculation
        public static BurnMode GetBurnMode(int index) { return burn_mode; }
        public static void SetBurnMode(BurnMode value) { burn_mode = value; }

        public static int GetLastManeuverIndex()
        {
            string vesselguid = GetVesselGuid();
            if (!HasFlightPlan())
                return -1;
            return plugin.FlightPlanNumberOfManoeuvres(vesselguid) - 1;
        }

        public static bool AddManeuver()
        {
            string vesselguid = GetVesselGuid();
            int last_burn_index = GetLastManeuverIndex();
            double initial_time = 0.0;
            bool inertially_fixed = false;

            EnsureFlightPlanExists();

            if (last_burn_index >= 0)
            {
                initial_time = GetManeuverTime(last_burn_index) + 60;
                inertially_fixed = GetManeuverIntertiallyFixed(last_burn_index);
            }
            else
            {
                initial_time = plugin.CurrentTime() + 60;
            }

            // TODO: fill in other parameters
            Burn candidate_burn = new Burn{
                thrust_in_kilonewtons = 10.0,
                specific_impulse_in_seconds_g0 = 100.0,
                frame = GenerateNavigationFrameParameters(),
                initial_time = initial_time,
                delta_v = new XYZ{x = 0.0,
                                  y = 0.0,
                                  z = 0.0},
                is_inertially_fixed = inertially_fixed};
            return plugin.FlightPlanAppend(vesselguid, candidate_burn);
        }

        public static void RemoveLastManeuver()
        {
            string vesselguid = GetVesselGuid();
            plugin.FlightPlanRemoveLast(vesselguid);
            if (GetLastManeuverIndex() < 0)
            {
                plugin.FlightPlanDelete(vesselguid);
            }
        }
    }
}  // namespace ksp_plugin_adapter
}  // namespace principia