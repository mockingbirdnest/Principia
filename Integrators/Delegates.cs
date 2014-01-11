using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Integrators {

  public delegate void ComputeRightHandSide(double[] y, double t, ref double[] result);

  public delegate void ComputeAutonomousRightHandSide(double[] y, ref double[] result);
}