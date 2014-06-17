namespace principia {
namespace rendering {

public delegate void RenderingFunction(float camera_distance);
public class MapRenderer : UnityEngine.MonoBehaviour {
  public static MapRenderer CreateAndAttach(RenderingFunction draw) {
    MapRenderer renderer
      = MapView.MapCamera.gameObject.GetComponent<MapRenderer>();
    if (renderer != null) {
      Destroy(renderer);
    }
    renderer = MapView.MapCamera.gameObject.AddComponent<MapRenderer>();
    renderer.draw_ = draw;
    return renderer;
  }

  private void OnPreCull() {
    draw_(MapView.MapCamera.Distance);
  }
  private RenderingFunction draw_;
}

}  // namespace rendering
}  // namespace principia
