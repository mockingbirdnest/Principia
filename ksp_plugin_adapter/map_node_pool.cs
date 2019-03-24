using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace principia {
namespace ksp_plugin_adapter {

internal class MapNodePool {

  public enum NodeSource {
    PREDICTION,
    FLIGHT_PLAN,
  }

  public MapNodePool() {
    nodes_ = new List<KSP.UI.Screens.Mapview.MapNode>();
    properties_ =
        new Dictionary<KSP.UI.Screens.Mapview.MapNode, MapNodeProperties>();
  }

  public void Clear() {
    foreach (var node in nodes_) {
      node.Terminate();
    }
    nodes_ = new List<KSP.UI.Screens.Mapview.MapNode>();
    properties_ =
        new Dictionary<KSP.UI.Screens.Mapview.MapNode, MapNodeProperties>();
    pool_index_ = 0;
  }

  public void Update() {
    for (int i = pool_index_; i < nodes_.Count; ++i) {
      if (properties_[nodes_[i]].visible) {
        properties_[nodes_[i]].visible = false;
        nodes_[i].NodeUpdate();
      }
    }
    for (int i = 0; i < pool_index_; ++i) {
      nodes_[i].NodeUpdate();
    }
    pool_index_ = 0;
  }

  public void RenderMarkers(DisposableIterator apsis_iterator,
                            MapObject.ObjectType type,
                            NodeSource source,
                            Vessel vessel,
                            CelestialBody celestial) {
    // We render at most 64 markers of one type and one provenance (e.g., at
    // most 64 perilunes for the prediction of the active vessel).  This is
    // more than is readable, and keeps the size of the map node pool under
    // control.
    for (int i = 0; i < 64 && !apsis_iterator.IteratorAtEnd();
         ++i, apsis_iterator.IteratorIncrement()) {
      QP apsis = apsis_iterator.IteratorGetDiscreteTrajectoryQP();
      MapNodeProperties node_properties = new MapNodeProperties {
        visible = true,
        object_type = type,
        vessel = vessel,
        celestial = celestial,
        world_position = (Vector3d)apsis.q,
        velocity = (Vector3d)apsis.p,
        source = source,
        time = apsis_iterator.IteratorGetDiscreteTrajectoryTime()
      };

      if (pool_index_ == nodes_.Count) {
        nodes_.Add(MakePoolNode());
      } else if (properties_[nodes_[pool_index_]].object_type != type) {
        // KSP attaches labels to its map nodes, but never detaches them.
        // If the node changes type, we end up with an arbitrary combination of
        // labels Ap, Pe, AN, DN.
        // Recreating the node entirely takes a long time (approximately 50 ns *
        // 𝑁, where 𝑁 is the total number of map nodes in existence), instead
        // we manually get rid of the labels.
        foreach (var component in
                 nodes_[pool_index_].transform.GetComponentsInChildren<
                     TMPro.TextMeshProUGUI>()) {
          UnityEngine.Object.Destroy(component.gameObject);
        }
      }
      properties_[nodes_[pool_index_++]] = node_properties;
    }
  }

  private KSP.UI.Screens.Mapview.MapNode MakePoolNode() {
    var new_node = KSP.UI.Screens.Mapview.MapNode.Create(
        "apsis",
        // If we see this colour, something has gone wrong.
        XKCDColors.Pale,
        pixelSize : 32,
        hoverable : true,
        pinnable : true,
        blocksInput : true);
    new_node.OnClick +=
        (KSP.UI.Screens.Mapview.MapNode node,
         Mouse.Buttons buttons) => {
          if (buttons == Mouse.Buttons.Left) {
            if (properties_[node].vessel != null) {
              if (PlanetariumCamera.fetch.target !=
                  properties_[node].vessel.mapObject) {
                PlanetariumCamera.fetch.SetTarget(
                    properties_[node].vessel.mapObject);
              }
            } else if (PlanetariumCamera.fetch.target !=
                       properties_[node].celestial.MapObject) {
              PlanetariumCamera.fetch.SetTarget(
                  properties_[node].celestial.MapObject);
            }
          }
        };
    new_node.OnUpdateVisible +=
        (KSP.UI.Screens.Mapview.MapNode node,
         KSP.UI.Screens.Mapview.MapNode.IconData icon) => {
          if (!properties_[node].visible) {
            icon.visible = false;
            return;
          }
          CelestialBody celestial = properties_[node].celestial;
          UnityEngine.Color colour =
              celestial.orbit == null
                  ? XKCDColors.SunshineYellow
                  : celestial.orbitDriver.Renderer.nodeColor;
          colour.a = 1;
          icon.visible = true;
          if (properties_[node].object_type == MapObject.ObjectType.Periapsis &&
              properties_[node].celestial.GetAltitude(
                  properties_[node].world_position) < 0) {
            // Make sure we see impacts.
            colour = XKCDColors.Orange;
          }
          if (properties_[node].object_type ==
              MapObject.ObjectType.ApproachIntersect) {
            colour = XKCDColors.Chartreuse;
          }
          icon.color = colour;
        };
    new_node.OnUpdateType +=
        (KSP.UI.Screens.Mapview.MapNode node,
         KSP.UI.Screens.Mapview.MapNode.TypeData type) => {
          if (properties_[node].object_type == MapObject.ObjectType.Periapsis &&
              properties_[node].celestial.GetAltitude(
                  properties_[node].world_position) < 0) {
            type.oType = MapObject.ObjectType.PatchTransition;
            type.pType =
                KSP.UI.Screens.Mapview.MapNode.PatchTransitionNodeType.Impact;
          } else if (properties_[node].object_type ==
                     MapObject.ObjectType.ApproachIntersect) {
            type.oType = properties_[node].object_type;
            type.aType = KSP.UI.Screens.Mapview.MapNode.ApproachNodeType
                             .CloseApproachOwn;
          } else {
            type.oType = properties_[node].object_type;
          }
        };
    new_node.OnUpdateCaption +=
        (KSP.UI.Screens.Mapview.MapNode node,
         KSP.UI.Screens.Mapview.MapNode.CaptionData caption) => {
          if (properties_[node].object_type == MapObject.ObjectType.Periapsis ||
              properties_[node].object_type == MapObject.ObjectType.Apoapsis) {
            String name =
                properties_[node].object_type == MapObject.ObjectType.Periapsis
                    ? "Periapsis"
                    : "Apoapsis";
            caption.Header =
                properties_[node].celestial.name + " " + name + " : <color=" +
                XKCDColors.HexFormat.Chartreuse + ">" +
                properties_[node].celestial.GetAltitude(
                    properties_[node].world_position).ToString(
                        "N0",
                        Culture.culture) + " m</color>";
            caption.captionLine2 =
                properties_[node].velocity.magnitude.ToString("N0",
                                                              Culture.culture) +
                " m/s";
          } else if (properties_[node].object_type ==
                         MapObject.ObjectType.AscendingNode ||
                     properties_[node].object_type ==
                         MapObject.ObjectType.DescendingNode) {
            String name = properties_[node].object_type ==
                                  MapObject.ObjectType.AscendingNode
                              ? "Ascending Node"
                              : "Descending Node";
            caption.Header = name;
            caption.captionLine2 =
                properties_[node].velocity.z.ToString("N0", Culture.culture) +
                " m/s out of plane";
          } else if (properties_[node].object_type ==
                     MapObject.ObjectType.ApproachIntersect) {
            caption.Header = "Target Approach : <color=" +
                             XKCDColors.HexFormat.Chartreuse + ">" +
                             (properties_[node].vessel.GetWorldPos3D() -
                              properties_[node].world_position)
                                 .magnitude.ToString("N0", Culture.culture) +
                             " m</color>";
            caption.captionLine2 =
                properties_[node].velocity.magnitude.ToString("N0",
                                                              Culture.culture) +
                " m/s";
          }
          caption.captionLine1 =
              "T" +
              FlightPlanner.FormatTimeSpan(TimeSpan.FromSeconds(
                  Planetarium.GetUniversalTime() - properties_[node].time));
          if (properties_[node].celestial.GetAltitude(
                  properties_[node].world_position) < 0 &&
              properties_[node].object_type == MapObject.ObjectType.Periapsis) {
            caption.Header = properties_[node].celestial.name +
                             " Impact<color=" +
                             XKCDColors.HexFormat.Chartreuse + "></color>";
            caption.captionLine1 = "";
            caption.captionLine2 = "";
          }
          switch (properties_[node].source) {
            case NodeSource.FLIGHT_PLAN:
              caption.Header = "Planned " + caption.Header;
              break;
            case NodeSource.PREDICTION:
              caption.Header = "Predicted " + caption.Header;
              break;
            default:
              throw Log.Fatal("Unexpected node source " +
                              properties_[node].source);
          }
        };
    new_node.OnUpdatePosition +=
        (KSP.UI.Screens.Mapview.MapNode node) =>
        ScaledSpace.LocalToScaledSpace(
            properties_[node].world_position);
    return new_node;
  }

  private class MapNodeProperties {
    public bool visible;
    public MapObject.ObjectType object_type;
    public Vector3d world_position;
     // Velocity in the plotting frame.  Note that the handedness is
     // inconsistent with World; for practical purposes only the norm or
     // individual coordinates of this vector should be used here.
    public Vector3d velocity;
    public Vessel vessel;
    public CelestialBody celestial;
    public NodeSource source;
    public double time;
  }

  private List<KSP.UI.Screens.Mapview.MapNode> nodes_;
  private Dictionary<KSP.UI.Screens.Mapview.MapNode,
                     MapNodeProperties> properties_;
  private int pool_index_ = 0;
}

}  // namespace ksp_plugin_adapter
}  // namespace principia
