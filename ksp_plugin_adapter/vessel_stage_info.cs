using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace principia { 
namespace ksp_plugin_adapter {
  public struct VesselStageInfo {
    public double start_mass;   // Tonnes
    public double end_mass;     // Tonnes
    public double thrust_vac;   // kN
    public double isp_vac;      // Seconds (g0)
    public double stage_deltav; // m/s

    public override string ToString()
    {
      return string.Format("start mass: {0} end mass: {1} thrust: {2} isp: {3} dv: {4}", start_mass, end_mass, thrust_vac, isp_vac, stage_deltav);
    }
  }  

  public interface IStageInfoProvider
  {
    string ProviderName();
    VesselStageInfo GetStageInfo(int stage_index);
  }

  public class StageInfoFactory
  {
    // Not implemented yet.
    static bool mj_installed_ = false;
    public static IStageInfoProvider Create(Vessel vessel) {
      if (mj_installed_){
        return new MechJebStageInfo(vessel);
      }

      //Mechjeb is not installed
      return new StockDeltaVAppStageInfo(vessel);
    }
  }

  // Not implemented yet.
  public class MechJebStageInfo : IStageInfoProvider
  {
    Vessel vessel_;
    public MechJebStageInfo(Vessel vessel){
      vessel_ = vessel;
    }

    public VesselStageInfo GetStageInfo(int stage_index)
    {
      throw new NotImplementedException();
    }

    public string ProviderName()
    {
      return "MechJeb 2";
    }
  }

  public class StockDeltaVAppStageInfo : IStageInfoProvider
  {
    private readonly Vessel vessel_;

    public StockDeltaVAppStageInfo(Vessel vessel){
      vessel_ = vessel;
    }

    public VesselStageInfo GetStageInfo(int stage_index) {
      if (stage_index < 0 || stage_index > vessel_.currentStage) {
        throw new ArgumentException("Invalid stage index");
      }

      DeltaVStageInfo stageInfo = vessel_.VesselDeltaV.GetStage(stage_index);
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
