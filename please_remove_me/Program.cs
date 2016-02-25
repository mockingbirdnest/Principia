using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.InteropServices;
using System.Text;
using System.Threading.Tasks;

namespace principia.ksp_plugin_adapter {
  class Program {
    static void Main(string[] args) {
      IntPtr plugin = Interface.NewPlugin(0.0, 0.0);
      int? i = null;
      plugin.InsertCelestialAbsoluteCartesian(0,
                                              i,
                                              "1m^3/s^2",
                                              "1m",
                                              "0 deg",
                                              "0 deg",
                                              null,
                                              null,
                                              "0 m",
                                              "0 m",
                                              "0 m",
                                              "0 m/s",
                                              "0m/s",
                                              "0m/s");
      i=0;
      plugin.InsertCelestialAbsoluteCartesian(1,
                                              i,
                                              "1m^3/s^2",
                                              "1m",
                                              "0 deg",
                                              "0 deg",
                                              null,
                                              null,
                                              "10 km",
                                              "0 m",
                                              "0 m",
                                              "0 m/s",
                                              "0m/s",
                                              "0m/s");
      i=1;
      plugin.InsertCelestialAbsoluteCartesian(2,
                                              i,
                                              "1m^3/s^2",
                                              "1m",
                                              "0 deg",
                                              "0 deg",
                                              null,
                                              null,
                                              "10 km",
                                              "0 m",
                                              "0 m",
                                              "0 m/s",
                                              "0m/s",
                                              "0m/s");
    }
  }
}
