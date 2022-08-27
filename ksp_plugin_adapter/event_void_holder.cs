using System.Collections.Generic;

namespace principia {
namespace ksp_plugin_adapter {

class EventVoidHolder {
  public void Add(EventVoid event_void, EventVoid.OnEvent on_event) {
    event_void.Add(on_event);
    registered_.Add(event_void, on_event);
  }

  public void RemoveAll() {
    foreach (var item in registered_) {
      item.Key.Remove(item.Value);
    }
    registered_.Clear();
  }

  private readonly Dictionary<EventVoid, EventVoid.OnEvent> registered_ =
      new Dictionary<EventVoid, EventVoid.OnEvent>();
}

} // namespace ksp_plugin_adapter
} // namespace principia
