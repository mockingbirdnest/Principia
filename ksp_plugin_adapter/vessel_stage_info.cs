using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

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

  public interface IStageInfoProvider {
    string ProviderName();
    public VesselStageInfo GetStageInfo(Vessel vessel, int stage_index);
  }

  public class StageInfoFactory {
    public static IStageInfoProvider Create() {
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
  public class MechJebStageInfo : IStageInfoProvider {
    public MechJebStageInfo(){
    }

    public VesselStageInfo GetStageInfo(Vessel vessel, int stage_index) {
      throw new NotImplementedException();
    }

    public string ProviderName() {
      return "MechJeb 2";
    }
  }

  public class StockDeltaVAppStageInfo : IStageInfoProvider {
    
    public VesselStageInfo GetStageInfo(Vessel vessel, int stage_index) {
      // Sadly, Stock Δv app won't update on RCS fuel being consumed,
      // and the update is asynchronous, so the user may need to
      // press the button twice. It's recommended to use MechJeb2 instead.
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

    public string ProviderName() {
      return "Stock DeltaV App";
    }
  }
} // namespace ksp_plugin_adapter
} // namespace principia
