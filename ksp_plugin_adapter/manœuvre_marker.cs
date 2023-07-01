using System;

namespace principia.ksp_plugin_adapter {

// Note that the 'base' scale of the marker is such that it is circumscribed
// by a sphere 1 unit in radius.
internal class ManœuvreMarker : UnityEngine.MonoBehaviour, IMouseEvents {
  public bool IsHovered { get; private set; } = false;
  public bool IsDragged { get; private set; } = false;
  public bool IsPinned { get; private set;} = false;
  public bool IsInteracting => IsHovered || IsDragged;

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

    transform.SetLayerRecursive((int)PrincipiaPluginAdapter.UnityLayers.
                                    Atmosphere);

    caption_ = MakeCaption();
    caption_.transform.SetParent(transform, false);
    caption_text_ = caption_.GetComponentInChildren<TMPro.TextMeshPro>();
    caption_.SetLayerRecursive((int)PrincipiaPluginAdapter.UnityLayers.UI);

    Disable();
  }

  public void Render(Vector3d world_position,
                     Vector3d initial_plotted_velocity,
                     int index,
                     NavigationManoeuvreFrenetTrihedron trihedron,
                     Func<NavigationManoeuvre> get_manœuvre,
                     Action<Burn> modify_burn) {
    initial_plotted_velocity_ = initial_plotted_velocity;
    index_ = index;
    get_manœuvre_ = get_manœuvre;
    modify_burn_ = modify_burn;

    transform.position = ScaledSpace.LocalToScaledSpace(world_position);

    tangent_.transform.localRotation = UnityEngine.Quaternion.FromToRotation(
        UnityEngine.Vector3.up,
        (Vector3d)trihedron.tangent);
    normal_.transform.localRotation = UnityEngine.Quaternion.FromToRotation(
        UnityEngine.Vector3.up,
        (Vector3d)trihedron.normal);
    binormal_.transform.localRotation = UnityEngine.Quaternion.FromToRotation(
        UnityEngine.Vector3.up,
        (Vector3d)trihedron.binormal);

    normalized_scale_ =
        (ScaledSpace.ScaledToLocalSpace(MapView.MapCamera.transform.position) -
         world_position).magnitude *
        0.015 *
        ScaledSpace.InverseScaleFactor;
    UpdateScale();

    UpdateColors();
    UpdateCaption();

    gameObject.SetActive(true);
  }

  public void Disable() {
    gameObject.SetActive(false);
    IsHovered = false;
    IsDragged = false;
  }

  private void SetColor(UnityEngine.GameObject go, UnityEngine.Color color) {
    if (IsInteracting) {
      UnityEngine.Color.RGBToHSV(color, out float h, out float s, out float v);
      color = UnityEngine.Color.HSVToRGB(h, s, v * 1.375f, hdr: true);
    }
    go.GetComponentInChildren<UnityEngine.Renderer>().material.color = color;
  }

  private void UpdateColors() {
    SetColor(base_, XKCDColors.LightGrey);
    SetColor(tangent_, Style.Tangent);
    SetColor(normal_, Style.Normal);
    SetColor(binormal_, Style.Binormal);
  }

  private void UpdateCaption() {
    if (!IsInteracting && !IsPinned) {
      caption_.SetActive(false);
      return;
    }
    caption_.SetActive(true);
    caption_.transform.position = ScaledToUIPosition(transform.position).position;

    if (IsInteracting) {
      caption_text_.color = caption_color;
    } else {
      caption_text_.color = caption_color_pinned;
    }

    var manœuvre = get_manœuvre_();
    var burn = manœuvre.burn;
    var caption = L10N.CacheFormat(
      "#Principia_MapNode_ManœuvreCaption",
      index_ + 1,
      ((Vector3d)burn.delta_v).magnitude.ToString("0.000"),
      manœuvre.duration.ToString("0.0"),
      "T" + new PrincipiaTimeSpan(
              Planetarium.GetUniversalTime() - burn.initial_time).
          Format(with_leading_zeroes: false, with_seconds: true).
          Replace(Culture.culture.NumberFormat.NegativeSign, "-")
    );

    caption_text_.SetText(caption);
  }

  private void UpdateScale() {
    hover_scale_offset_ = UnityEngine.Mathf.MoveTowards(hover_scale_offset_,
      IsInteracting ? 1.25f : 1f,
      0.075f);
    var resulting_scale = normalized_scale_ * hover_scale_offset_;
    transform.localScale = Vector3d.one * resulting_scale;
    // Inversely scale the label, since it is drawn in UI coordinates and thus
    // does not require additional dynamic scaling.
    caption_.transform.localScale = Vector3d.one / resulting_scale;
  }

  private UnityEngine.Vector3 ScreenManœuvrePosition() =>
      ScaledToFlattenedScreenPosition(transform.position);

  public void OnMouseEnter() {
    if (!isActiveAndEnabled) {
      return;
    }
    IsHovered = true;
  }

  public void OnMouseDown() {
    mouse_offset_at_click_ =
        UnityEngine.Input.mousePosition - ScreenManœuvrePosition();
    IsDragged = true;
  }

  public void OnMouseOver() {
    if (!isActiveAndEnabled) {
      return;
    }
    if (Mouse.Right.GetButtonUp()) {
      IsPinned = !IsPinned;
    }
  }

  public void OnMouseDrag() {
    if (!isActiveAndEnabled) {
      return;
    }

    var mouse_offset_now =
        UnityEngine.Input.mousePosition - ScreenManœuvrePosition();
    var mouse_displacement = mouse_offset_now - mouse_offset_at_click_;

    // TODO(egg): This is cheesy. Could we do it in the C++?
    var screen_velocity =
        (ScaledToFlattenedScreenPosition(
             transform.position +
             ScaledSpace.InverseScaleFactor * initial_plotted_velocity_) -
         ScreenManœuvrePosition());

    var Δt = UnityEngine.Vector3.Dot(screen_velocity, mouse_displacement) /
              screen_velocity.sqrMagnitude;
    var burn = get_manœuvre_().burn;
    burn.initial_time += Δt;
    modify_burn_(burn);
  }

  public void OnMouseUp() {
    if (!isActiveAndEnabled) {
      return;
    }
    IsDragged = false;
  }

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

  private static UnityEngine.GameObject MakeCaption() {
    var caption = new UnityEngine.GameObject("caption");
    var caption_inner = new UnityEngine.GameObject("TMPro");
    var text = caption_inner.AddComponent<TMPro.TextMeshPro>();
    var inner_transform = caption_inner.GetComponent<UnityEngine.RectTransform>();
    inner_transform.SetParent(caption.transform, false);
    inner_transform.localPosition = new UnityEngine.Vector3(0f, -10f, 0f);
    inner_transform.sizeDelta = new UnityEngine.Vector2(300f, 1f);
    text.font = UISkinManager.TMPFont;
    text.fontSize = 140;
    text.alignment = TMPro.TextAlignmentOptions.Top;
    return caption;
  }

  private static UnityEngine.Vector3 ScaledToFlattenedScreenPosition(
      Vector3d scaled_position) {
    var position = PlanetariumCamera.Camera.WorldToScreenPoint(scaled_position);
    position.z = 0;
    return position;
  }

  private static (UnityEngine.Vector3 position, bool visible) ScaledToUIPosition(
    Vector3d scaled_position) {
      bool visible = false;
      var position = KSP.UI.Screens.Mapview.MapViewCanvasUtil.ScaledToUISpacePos(
          scaled_position,
          ref visible,
          KSP.UI.Screens.Mapview.MapNode.zSpaceEasing,
          KSP.UI.Screens.Mapview.MapNode.zSpaceMidpoint,
          KSP.UI.Screens.Mapview.MapNode.zSpaceUIStart,
          KSP.UI.Screens.Mapview.MapNode.zSpaceLength);
    return (position, visible);
  }

  private UnityEngine.GameObject base_;
  private UnityEngine.GameObject tangent_;
  private UnityEngine.GameObject normal_;
  private UnityEngine.GameObject binormal_;

  private UnityEngine.GameObject caption_;
  private TMPro.TextMeshPro caption_text_;

  private double normalized_scale_;
  private float hover_scale_offset_;

  private UnityEngine.Vector3 mouse_offset_at_click_;

  private Vector3d initial_plotted_velocity_;
  private int index_;
  private Func<NavigationManoeuvre> get_manœuvre_;
  private Action<Burn> modify_burn_;

  private static UnityEngine.Material marker_material_;

  public static UnityEngine.Material marker_material {
    get {
      if (marker_material_ == null) {
        marker_material_ =
            new UnityEngine.Material(UnityEngine.Shader.Find("Unlit/Color"));
      }
      return marker_material_;
    }
  }

  public static UnityEngine.Color caption_color
    = new UnityEngine.Color(191f / 255f, 1f, 0f, 1f);
  public static UnityEngine.Color caption_color_pinned
    = new UnityEngine.Color(191f / 255f, 1f, 0f, 0.6f);
}

}
