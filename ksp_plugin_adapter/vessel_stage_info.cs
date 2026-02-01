using System;
using System.Collections.Generic;
using System.Linq;

namespace principia { 
namespace ksp_plugin_adapter {
  public class VesselStageInfo {
    public double start_mass;   // Tonnes
    public double end_mass;     // Tonnes
    public double thrust_vac;   // kN
    public double isp_vac;      // Seconds (g0)
    public double stage_deltav; // m/s

    public override string ToString() {
      return string.Format("start mass: {0} end mass: {1} thrust: {2} isp: {3} dv: {4}", start_mass, end_mass, thrust_vac, isp_vac, stage_deltav);
    }
  }

  public abstract class StageInfoProviderBase {
    public abstract string ProviderName { get; }
    public abstract VesselStageInfo GetStageInfo(Vessel vessel, int stage_index);

    public VesselStageInfo GetActiveEngines(Vessel vessel) {
      ModuleEngines[] active_engines =
        (from part in vessel.parts
         select (from PartModule module in part.Modules
                 where (module as ModuleEngines)?.EngineIgnited == true
                 select module as ModuleEngines)).SelectMany(x => x).ToArray();
      Vector3d reference_direction = vessel.ReferenceTransform.up;
      double[] thrusts =
          (from engine in active_engines
           select
               engine.MaxThrustOutputVac(useThrustLimiter: true) *
               (from transform in engine.thrustTransforms
                select Math.Max(0,
                                Vector3d.Dot(reference_direction,
                                             -transform.forward))).Average()).
          ToArray();
      double thrust_in_kilonewtons_ = thrusts.Sum();

      // This would use zip if we had 4.0 or later.  We loop for now.
      double Σ_f_over_i_sp = 
                 thrusts.Zip(active_engines, 
                             (thrust, engine) 
                           => thrust / engine.atmosphereCurve
                                             .Evaluate(0)).Sum();

      /*for (int i = 0; i < active_engines.Length; ++i) {
        Σ_f_over_i_sp +=
            thrusts[i] / active_engines[i].atmosphereCurve.Evaluate(0);
      }*/

      double specific_impulse_in_seconds_g0_ = thrust_in_kilonewtons_ / Σ_f_over_i_sp;

      return new VesselStageInfo {
        start_mass = vessel.totalMass, // Not used
        end_mass = 0, // Not used
        thrust_vac = thrust_in_kilonewtons_,
        isp_vac = specific_impulse_in_seconds_g0_,
        stage_deltav = double.PositiveInfinity // Not used
      };
    }

    public VesselStageInfo GetActiveRCS(Vessel vessel) {
      ModuleRCS[] active_rcs = (from part in vessel.parts
                              select (from PartModule module in part.Modules
                                      where module is ModuleRCS module_rcs &&
                                            module_rcs.rcsEnabled
                                      select module as ModuleRCS)).
        SelectMany(x => x).ToArray();
      Vector3d reference_direction = vessel.ReferenceTransform.up;
      // NOTE(egg): NathanKell informs me that in >= 1.0.5, RCS has a useZaxis
      // property, that controls whether they thrust -up or -forward.  The madness
      // keeps piling up.
      double[] thrusts = (from engine in active_rcs
                          select engine.thrusterPower *
                                 (from transform in engine.thrusterTransforms
                                  where transform.gameObject.activeInHierarchy
                                  select Math.Max(0,
                                                  Vector3d.Dot(
                                                      reference_direction,
                                                      -transform.up))).Sum()).
          ToArray();
      double thrust_in_kilonewtons_ = thrusts.Sum();

      // This would use zip if we had 4.0 or later.  We loop for now.
      double Σ_f_over_i_sp = 
                 thrusts.Zip(active_rcs,
                             (thrust, engine)
                           => thrust / engine.atmosphereCurve
                                             .Evaluate(0)).Sum();
      /*
      for (int i = 0; i < active_rcs.Length; ++i) {
        Σ_f_over_i_sp +=
            thrusts[i] / active_rcs[i].atmosphereCurve.Evaluate(0);
      }*/
      double specific_impulse_in_seconds_g0_ = thrust_in_kilonewtons_ / Σ_f_over_i_sp;

      return new VesselStageInfo {
        start_mass = vessel.totalMass, // Not used
        end_mass = 0, // Not used
        thrust_vac = thrust_in_kilonewtons_,
        isp_vac = specific_impulse_in_seconds_g0_,
        stage_deltav = double.PositiveInfinity // Not used
      };
    }
  }

  public class StageInfoFactory {
    public static StageInfoProviderBase Create() {
      // Prioritize stage info provided by MechJeb
      if (mechjeb_installed()) {
        return new MechJebStageInfo();
      }

      // Mechjeb is not installed
      return new StockDeltaVAppStageInfo();
    }

    // Not implemented yet.
    private static bool mechjeb_installed() {
      return false;
    }
  }

  // Not implemented yet.
  public class MechJebStageInfo : StageInfoProviderBase {
    public MechJebStageInfo(){
    }

    public override VesselStageInfo GetStageInfo(Vessel vessel, int stage_index) {
      throw new NotImplementedException();
    }

    public override string ProviderName => "MechJeb 2";
  }

  public class StockDeltaVAppStageInfo : StageInfoProviderBase {
    
    public override VesselStageInfo GetStageInfo(Vessel vessel, int stage_index) {
      // Sadly, Stock Δv app won't update on RCS fuel being consumed,
      // and the update is asynchronous, so the user may need to
      // press the button twice to get the correct stage mass.
      // It's recommended to use MechJeb2 instead.
      vessel.VesselDeltaV.SetCalcsDirty(resetPartCaches: false);
      DeltaVStageInfo stageInfo = vessel.VesselDeltaV.GetStage(stage_index);

      if (stageInfo == null) {
        return null;
      }

      return new VesselStageInfo{
        start_mass = stageInfo.startMass,
        end_mass = stageInfo.endMass,
        thrust_vac = stageInfo.thrustVac,
        isp_vac = stageInfo.ispVac,
        stage_deltav = stageInfo.deltaVinVac
      };
    }

    public override string ProviderName => "Stock DeltaV App";
  }
} // namespace ksp_plugin_adapter
} // namespace principia
