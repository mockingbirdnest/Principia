using System;
using System.Collections.Generic;
using System.Security.Authentication;

namespace principia {
namespace ksp_plugin_adapter {

public class ExternalInterface {
  public XYZ CelestialGetPosition(int body_index, double time) {
    ThrowOnError(
        adapter_.Plugin().ExternalCelestialGetPosition(
            body_index, time, out XYZ result));
    return result;
  }

  public XYZ CelestialGetSurfacePosition(
      int body_index,
      double planetocentric_latitude_in_degrees,
      double planetocentric_longitude_in_degrees,
      double radius,
      double time) {
    ThrowOnError(
        adapter_.Plugin().ExternalCelestialGetSurfacePosition(
            body_index, planetocentric_latitude_in_degrees,
            planetocentric_longitude_in_degrees, radius, time, out XYZ result));
    return result;
  }

  public XY GeopotentialGetCoefficient(int body_index, int degree, int order) {
    ThrowOnError(
        adapter_.Plugin().ExternalGeopotentialGetCoefficient(
            body_index, degree, order, out XY result));
    return result;
  }

  public double GeopotentialGetReferenceRadius(int body_index) {
    ThrowOnError(
        adapter_.Plugin().ExternalGeopotentialGetReferenceRadius(
            body_index, out double result));
    return result;
  }

  public XYZ VesselGetPosition(string vessel_guid, double time) {
    ThrowOnError(
        adapter_.Plugin().ExternalVesselGetPosition(
            vessel_guid, time, out XYZ result));
    return result;
  }

  public static ExternalInterface Get() {
    List<ScenarioModule> modules;
    try {
      modules = ScenarioRunner.GetLoadedModules();
    } catch {
      return null;
    }
    foreach (var module in modules) {
      if (module is PrincipiaPluginAdapter adapter) {
        return new ExternalInterface(adapter);
      }
    }
    return null;
  }

  private static void ThrowOnError(Status status) {
    switch (status.error) {
      case 0:
        return;
      case 1:
        throw new OperationCanceledException($"CANCELLED: {status.message}");
      case 2:
        throw new Exception($"UNKNOWN: {status.message}");
      case 3:
        throw new ArgumentException($"INVALID_ARGUMENT: {status.message}");
      case 4:
        throw new TimeoutException($"DEADLINE_EXCEEDED: {status.message}");
      case 5:
        throw new KeyNotFoundException($"NOT_FOUND: {status.message}");
      case 6:
        throw new ArgumentException($"ALREADY_EXISTS: {status.message}");
      case 7:
        throw new UnauthorizedAccessException(
            $"PERMISSION_DENIED: {status.message}");
      case 16:
        throw new AuthenticationException($"UNAUTHENTICATED: {status.message}");
      case 8:
        throw new Exception($"RESOURCE_EXHAUSTED: {status.message}");
      case 9:
        throw new Exception($"FAILED_PRECONDITION: {status.message}");
      case 10:
        throw new Exception($"ABORTED: {status.message}");
      case 11:
        throw new ArgumentOutOfRangeException(
            $"OUT_OF_RANGE: {status.message}");
      case 12:
        throw new NotImplementedException($"UNIMPLEMENTED: {status.message}");
      case 13:
        throw new Exception($"INTERNAL: {status.message}");
      case 14:
        throw new Exception($"UNAVAILABLE: {status.message}");
      case 15:
        throw new Exception($"DATA_LOSS: {status.message}");
      default:
        throw new Exception($"Error {status.error}: {status.message}");
    }
  }

  private ExternalInterface(PrincipiaPluginAdapter adapter) {
    adapter_ = adapter;
    if (!Loader.loaded_principia_dll_) {
      throw new DllNotFoundException(
          "The Principia native DLL failed to load");
    }
  }

  private readonly PrincipiaPluginAdapter adapter_;
}

}  // namespace ksp_plugin_adapter
}  // namespace principia
