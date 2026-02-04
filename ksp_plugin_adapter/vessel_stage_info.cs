using System;
using System.Collections;
using System.Linq;
using System.Reflection;
using UnityEngine;

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
    public abstract bool IsAvailable { get; }
    public abstract string ProviderName { get; }
    public abstract VesselStageInfo GetStageInfo(Vessel vessel, int stage_index);

    public virtual VesselStageInfo GetActiveEngines(Vessel vessel) {
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

    public virtual VesselStageInfo GetActiveRCS(Vessel vessel) {
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
      var mechjeb_provider = new MechJebStageInfo();
      if (mechjeb_provider.IsAvailable) {
        return mechjeb_provider;
      }

      // Mechjeb is not installed, fallback to StockDeltaVApp
      return new StockDeltaVAppStageInfo();
    }
  }

  public class MechJebStageInfo : StageInfoProviderBase {
    public override bool IsAvailable => is_available_;
    public override string ProviderName => "MechJeb 2";

    private bool is_available_;

    private readonly Type MechJebCore_;
    private readonly Type ModuleStageStats_;
    private readonly Type FuelStats_;
    private readonly Type VesselExtensions_;

    private readonly MethodInfo GetMasterMechJeb_;
    private readonly MethodInfo GetComputerModule_;

    private readonly FieldInfo VacStats_;
    private readonly FieldInfo StartMass_;
    private readonly FieldInfo EndMass_;
    private readonly FieldInfo Thrust_;
    private readonly FieldInfo Isp_;
    private readonly FieldInfo DeltaV_;

    public MechJebStageInfo(){
      try {
        MechJebCore_ = Type.GetType("MuMech.MechJebCore, MechJeb2", true);
        VesselExtensions_ = Type.GetType("MuMech.VesselExtensions, MechJeb2", true);
        ModuleStageStats_ = Type.GetType("MuMech.MechJebModuleStageStats, MechJeb2", true);
        FuelStats_ = Type.GetType("MechJebLib.FuelFlowSimulation.FuelStats, MechJebLib", true);

        GetMasterMechJeb_ = GetMethod_(VesselExtensions_, "GetMasterMechJeb", [typeof(Vessel)]);
        GetComputerModule_ = GetMethod_(MechJebCore_, "GetComputerModule", [])
                            .MakeGenericMethod(ModuleStageStats_);

        VacStats_  = GetField_(ModuleStageStats_, "VacStats");

        StartMass_ = GetField_(FuelStats_, "StartMass");
        EndMass_   = GetField_(FuelStats_, "EndMass");
        Thrust_    = GetField_(FuelStats_, "Thrust");
        Isp_       = GetField_(FuelStats_, "Isp");
        DeltaV_    = GetField_(FuelStats_, "DeltaV");

        is_available_ = true;
      } catch (Exception e) {
        is_available_ = false;
        Debug.LogWarning("[Principia] Failed to load MechJeb2 Methods");
        Debug.LogException(e);
      }
    }

    public override VesselStageInfo GetStageInfo(Vessel vessel, int stage_index) {
      if (!IsAvailable) {
        return null;
      }

      var mechjeb_ = GetMasterMechJeb_.Invoke(null, [vessel]);
      if (mechjeb_ == null) {
        Debug.LogWarning("[Principia] Cannot get MechJeb Core");
        return null;
      }

      var stage_stats_module = GetComputerModule_.Invoke(mechjeb_, []);
      if (stage_stats_module == null) {
        Debug.LogWarning("[Principia] Cannot get MechJeb StageStats module");
        return null;
      }

      var stages = (IList) VacStats_.GetValue(stage_stats_module);

      if (stage_index < 0 || stage_index >= stages.Count) {
        Debug.LogWarning("[Principia] Invalid stage: " + stage_index.ToString());
        return null;
      }

      var stage_info = stages[stage_index];

      return new VesselStageInfo {
        start_mass   = (double) StartMass_.GetValue(stage_info),
        end_mass     = (double) EndMass_.GetValue(stage_info),
        thrust_vac   = (double) Thrust_.GetValue(stage_info),
        isp_vac      = (double) Isp_.GetValue(stage_info),
        stage_deltav = (double) DeltaV_.GetValue(stage_info)
      };
    }

    private MethodInfo GetMethod_(Type type, string name, Type[] types) {
      return type.GetMethod(name, types) ?? throw new MissingMethodException(type.Name, name);
    }

    private FieldInfo GetField_(Type type, string name) {
      return type.GetField(name) ?? throw new MissingFieldException(type.Name, name);
    }
  }

  public class StockDeltaVAppStageInfo : StageInfoProviderBase {
    public override bool IsAvailable => true;
    public override string ProviderName => "Stock DeltaV App";
    
    public override VesselStageInfo GetStageInfo(Vessel vessel, int stage_index) {
      // Sadly, Stock Δv app won't update on RCS fuel being consumed,
      // and the update is asynchronous, so the user may need to
      // press the button twice to get the correct stage mass.
      // It's recommended to use MechJeb2 instead.
      vessel.VesselDeltaV.SetCalcsDirty(resetPartCaches: false);
      var stage_info = vessel.VesselDeltaV.GetStage(stage_index);

      if (stage_info == null) {
        Debug.LogWarning("[Principia] Invalid stage: " + stage_index.ToString());
        return null;
      }

      return new VesselStageInfo{
        start_mass = stage_info.startMass,
        end_mass = stage_info.endMass,
        thrust_vac = stage_info.thrustVac,
        isp_vac = stage_info.ispVac,
        stage_deltav = stage_info.deltaVinVac
      };
    }
  }
} // namespace ksp_plugin_adapter
} // namespace principia
