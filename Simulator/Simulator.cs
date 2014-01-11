using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Simulator {

  internal class Simulator {

    private static void Main(string[] args) {
      SimulationParameters parameters = new SimulationParameters(args[0]);
      NBodySystem system = new NBodySystem(parameters.m, parameters.q0, parameters.v0);
      Export(system.Simulate(parameters.tmax, parameters.Δt, parameters.samplingPeriod), args[1]);
    }

    private static void Export(Integrators.SymplecticPartitionedRungeKutta.Solution simulation, string path) {
      int N = simulation.Momentum[0].Length / 3;
      StreamWriter file = new StreamWriter(path);
      file.WriteLine("{{"); // Position.
      for (int i = 0; i < N; ++i) {
        file.WriteLine("{"); // Body i.
        for (int t = 0; t < simulation.Time.Length; ++t) {
          file.Write("{"); // Timestep t.
          file.Write(simulation.Position[t][3 * i] + ", "
            + simulation.Position[t][3 * i + 1] + ", "
            + simulation.Position[t][3 * i + 2]);
          file.Write("}"); // End timestep t.
          if (t + 1 < simulation.Time.Length) {
            file.WriteLine(",");
          }
        }
        file.WriteLine("}"); // End body i.
        if (i + 1 < N) {
          file.Write(",");
        }
      }
      file.WriteLine("}, {"); // Velocity.
      for (int i = 0; i < N; ++i) {
        file.WriteLine("{"); // Body i.
        for (int t = 0; t < simulation.Time.Length; ++t) {
          file.Write("{"); // Timestep t.
          file.Write(simulation.Momentum[t][3 * i] + ", "
            + simulation.Momentum[t][3 * i + 1] + ", "
            + simulation.Momentum[t][3 * i + 2]);
          file.Write("}"); // End timestep t.
          if (t + 1 < simulation.Time.Length) {
            file.WriteLine(",");
          }
        }
        file.WriteLine("}"); // End body i.
        if (i + 1 < N) {
          file.Write(",");
        }
      }
      file.WriteLine("}, {"); // Time.
      for (int t = 0; t < simulation.Time.Length; ++t) {
        file.Write(simulation.Time[t]);
        if (t + 1 < simulation.Time.Length) {
          file.WriteLine(",");
        }
      }
      file.WriteLine("}}"); // The end.
      file.Close();
    }

    private class SimulationParameters {

      // Expected file format (one number per line):
      // - number of bodies N (int)
      // - masses (double[N])
      // - q0 (double[3N])
      // - v0 (double[3N])
      // - tmax (double)
      // - Δt (double)
      // - samplingPeriod (int)
      public SimulationParameters(string path) {
        StreamReader file = new StreamReader(path);
        int N = Int32.Parse(file.ReadLine());
        m = new double[N];
        for (int i = 0; i < N; ++i) {
          m[i] = Double.Parse(file.ReadLine());
        }
        q0 = new double[3 * N];
        for (int i = 0; i < 3 * N; ++i) {
          q0[i] = Double.Parse(file.ReadLine());
        }
        v0 = new double[3 * N];
        for (int i = 0; i < 3 * N; ++i) {
          v0[i] = Double.Parse(file.ReadLine());
        }
        tmax = Double.Parse(file.ReadLine());
        Δt = Double.Parse(file.ReadLine());
        samplingPeriod = Int32.Parse(file.ReadLine());
      }

      public double[] m;
      public double[] q0;
      public double[] v0;
      public double tmax;
      public double Δt;
      public int samplingPeriod;
    }
  }
}