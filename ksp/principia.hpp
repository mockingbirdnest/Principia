
namespace principia {
namespace ksp {

public ref class Principia : public UnityEngine::MonoBehaviour {
 private:
  // Unity event functions, accessed through reflection.
  // See http://docs.unity3d.com/Manual/ExecutionOrder.html.
  // Before the first frame update.
  void Start();
  // Updates.
  void FixedUpdate();
  void Update();
  // Rendering.
  void OnPreCull();
  // GUI rendering.
  void OnGUI();
  // When the object is destroyed.
  void OnDestroy();
};

}  // namespace ksp
}  // namespace principia
