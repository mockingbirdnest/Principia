using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using KSP.Localization;

namespace principia {
namespace ksp_plugin_adapter {

[KSPScenario(createOptions: ScenarioCreationOptions.AddToAllGames,
             tgtScenes: new GameScenes[]{GameScenes.FLIGHT,
                                         GameScenes.MAINMENU,
                                         GameScenes.SPACECENTER,
                                         GameScenes.TRACKSTATION})]
public partial class PrincipiaPluginAdapterStub : ScenarioModule,
                                                  SupervisedWindowRenderer.
                                                  ISupervisor {
  [KSPField(isPersistant = true)]
  private readonly Dialog dll_stub_executed_ = new Dialog(persist_state: false);

#pragma warning disable 67
  public event Action LockClearing;
  public event Action WindowsDisposal;
  public event Action WindowsRendering;
#pragma warning restore 67

  PrincipiaPluginAdapterStub() {
    dll_stub_executed_.message =
        Localizer.Format("#Principia_DLLStubExecuted");
    dll_stub_executed_.Show();
  }

  private void OnGUI() {
    dll_stub_executed_.RenderWindow();
  }
}

}  // namespace ksp_plugin_adapter
}  // namespace principia
