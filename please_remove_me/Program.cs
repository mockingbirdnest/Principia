using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.InteropServices;
using System.Text;
using System.Threading.Tasks;

namespace principia.ksp_plugin_adapter {
  class Program {
    static void Main(string[] args) {
      string s = Interface.SayHello();
      Console.WriteLine(s);
      Log.Error(Interface.SayHello());
      string version, build_date;
      Interface.GetVersion(out build_date, out version);
    }
  }
}
