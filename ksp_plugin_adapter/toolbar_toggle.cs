using System;
using System.Collections.Generic;

namespace principia {
namespace ksp_plugin_adapter {
  
  /// <summary>
  /// Represents a Blizzy's toolbar button that is toggled on/off.
  /// </summary>
  internal class ToolbarToggle {
    public ToolbarToggle(
        string id,
        string tooltip,
        GameScenesVisibility visibility,
        Func<bool> get = null,
        Action<bool> set = null,
        bool initialValue = false) {
      Id = id;
      Tooltip = tooltip;
      Visibility = visibility;
      get_ = get ?? (() => value_);
      set_ = set ?? (value => { });
      value_ = initialValue;
      
      if (toggles_.ContainsKey(id)) {
        toggles_[id].Desstroy();
        toggles_.Remove(id);
      }
      toggles_.Add(id, this);
      
      if (is_toolbar_available_) {
        toolbar_button_ = manager_.add(nameof(principia), nameof(principia) + id);
        toolbar_button_.ToolTip = tooltip;
        toolbar_button_.Visibility = visibility;
        toolbar_button_.OnClick += OnClick;
        SetTexture(Value);
      }
    }
    
    public string Id { get; }
    public string Tooltip { get; }
    public GameScenesVisibility Visibility { get; }
    
    public bool Value {
      get => get_();
      set {
        value_ = value;
        set_(value_);
        SetTexture(value_);
      }
    }

    public bool IsEnabled
    {
      get => is_enabled_;
      set {
        is_enabled_ = value;
        if (toolbar_button_ != null) {
          toolbar_button_.Enabled = is_enabled_;
        }
      }
    }

    public void Desstroy() {
      if (toolbar_button_ != null) {
        toolbar_button_.OnClick -= OnClick;
        toolbar_button_.Destroy();
      }
    }
    
    private bool value_;
    private bool is_enabled_;
    
    private readonly Func<bool> get_;
    private readonly Action<bool> set_;
    private readonly IButton toolbar_button_;
    
    private void SetTexture(bool value) {
      if (toolbar_button_ != null) {
        toolbar_button_.TexturePath = button_texture_path_ + Id + (value ? "-on" : "-off");
      }
    }
    
    private void OnClick(ClickEvent clickEvent) {
      Value = !Value;
    }
    
    private static readonly IDictionary<string, ToolbarToggle> toggles_ = new Dictionary<string, ToolbarToggle>();
    
    private static readonly bool is_toolbar_available_ = ToolbarManager.ToolbarAvailable;
    private static readonly IToolbarManager manager_ = ToolbarManager.Instance;
    private const string button_texture_path_ = "Principia/assets/toolbar/";
  }
  
}
}