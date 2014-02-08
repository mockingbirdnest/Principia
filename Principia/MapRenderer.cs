using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using UnityEngine;

namespace Principia {
  public delegate void RenderingFunction();
  public class MapRenderer : MonoBehaviour {
    private RenderingFunction draw;
    public static MapRenderer CreateAndAttach(RenderingFunction draw) {
      MapRenderer renderer
        = MapView.MapCamera.gameObject.GetComponent<MapRenderer>();
      if (renderer) {
        Destroy(renderer);
      }
      renderer = MapView.MapCamera.gameObject.AddComponent<MapRenderer>();
      renderer.draw = draw;
      return renderer;
    }
    public void OnPreCull() {
      draw();
    }
  }
}