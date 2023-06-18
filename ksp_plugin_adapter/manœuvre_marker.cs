namespace principia.ksp_plugin_adapter
{
  // Note that the 'base' scale of the marker is such that it is circumscribed
  // by a sphere 1 unit in radius.
  internal class ManœuvreMarker : UnityEngine.MonoBehaviour,
                                  IMouseEvents {
    public bool IsHovered {
      get;
      private set;
    } = false;

    public void Awake() {
      base_ = MakeTrihedronBase();
      tangent_ = MakeTrihedronArrow();
      normal_ = MakeTrihedronArrow();
      binormal_ = MakeTrihedronArrow();

      base_.transform.SetParent(gameObject.transform, false);
      tangent_.transform.SetParent(gameObject.transform, false);
      normal_.transform.SetParent(gameObject.transform, false);
      binormal_.transform.SetParent(gameObject.transform, false);

      var collider = gameObject.AddComponent<UnityEngine.SphereCollider>();
      collider.radius = 0.75f;

      transform.SetLayerRecursive((int)PrincipiaPluginAdapter.UnityLayers.Atmosphere);

      Disable();
    }

    public void Render(Vector3d world_position,
                       NavigationManoeuvreFrenetTrihedron manœuvre,
                       double scale) {
      transform.position = ScaledSpace.LocalToScaledSpace(world_position);

      tangent_.transform.localRotation = UnityEngine.Quaternion.FromToRotation(
          UnityEngine.Vector3.up,
          (Vector3d)manœuvre.tangent);
      normal_.transform.localRotation = UnityEngine.Quaternion.FromToRotation(
          UnityEngine.Vector3.up,
          (Vector3d)manœuvre.normal);
      binormal_.transform.localRotation = UnityEngine.Quaternion.FromToRotation(
          UnityEngine.Vector3.up,
          (Vector3d)manœuvre.binormal);

      transform.localScale = UnityEngine.Vector3.one * (float)scale
          * ScaledSpace.InverseScaleFactor;

      UpdateMarkerColors();

      gameObject.SetActive(true);
    }

    public void Disable() {
      gameObject.SetActive(false);
      IsHovered = false;
    }

    private void SetColor(UnityEngine.GameObject go, UnityEngine.Color color) {
      if (IsHovered) {
        UnityEngine.Color.RGBToHSV(color, out float h, out float s, out float v);
        color = UnityEngine.Color.HSVToRGB(h, s, v * 1.375f, hdr: true);
      }
      go.GetComponentInChildren<UnityEngine.Renderer>().material.color = color;
    }

    private void UpdateMarkerColors() {
      SetColor(base_, XKCDColors.LightGrey);
      SetColor(tangent_, Style.Tangent);
      SetColor(normal_, Style.Normal);
      SetColor(binormal_, Style.Binormal);
    }

    public void OnMouseEnter() {
      if (!isActiveAndEnabled) {
        return;
      }
      IsHovered = true;
    }

    public void OnMouseDown() {}

    public void OnMouseDrag() {}

    public void OnMouseUp() {}

    public void OnMouseExit() {
      if (!isActiveAndEnabled) {
        return;
      }
      IsHovered = false;
    }

    public UnityEngine.MonoBehaviour GetInstance() {
      return this;
    }

    private static UnityEngine.GameObject MakeTrihedronBase() {
      var sphere = UnityEngine.GameObject.CreatePrimitive(
          UnityEngine.PrimitiveType.Sphere);
      Destroy(sphere.GetComponent<UnityEngine.Collider>());
      sphere.transform.localScale *= 0.3f;
      sphere.GetComponent<UnityEngine.Renderer>().material =
        new UnityEngine.Material(marker_material);
      return sphere;
    }

    private static UnityEngine.GameObject MakeTrihedronArrow() {
      var arrow = new UnityEngine.GameObject("arrow");
      var cylinder = UnityEngine.GameObject.CreatePrimitive(
            UnityEngine.PrimitiveType.Cylinder);
      Destroy(cylinder.GetComponent<UnityEngine.Collider>());
      cylinder.transform.SetParent(arrow.transform, false);
      cylinder.transform.localScale = new UnityEngine.Vector3(0.15f, 0.5f, 0.15f);
      cylinder.transform.localPosition = new UnityEngine.Vector3(0f, 0.5f, 0f);
      cylinder.GetComponent<UnityEngine.Renderer>().material =
        new UnityEngine.Material(marker_material);
      return arrow;
    }

    private UnityEngine.GameObject base_;
    private UnityEngine.GameObject tangent_;
    private UnityEngine.GameObject normal_;
    private UnityEngine.GameObject binormal_;

    private static UnityEngine.Material marker_material_;
    public static UnityEngine.Material marker_material {
      get {
        if (marker_material_ == null) {
          marker_material_ = new UnityEngine.Material(UnityEngine.Shader.Find("Unlit/Color"));
        }
        return marker_material_;
        }
    }
  }
}
