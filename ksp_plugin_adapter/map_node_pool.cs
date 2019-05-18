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
                            ReferenceFrameSelector reference_frame) {
    // We render at most 64 markers of one type and one provenance (e.g., at
    // most 64 perilunes for the prediction of the active vessel).  This is
    // more than is readable, and keeps the size of the map node pool under
    // control.
    MapObject associated_map_object;
    UnityEngine.Color colour;
    if (type == MapObject.ObjectType.Apoapsis ||
        type == MapObject.ObjectType.Periapsis) {
      CelestialBody fixed_body = reference_frame.selected_celestial;
      associated_map_object = fixed_body.MapObject;
      colour = fixed_body.orbit == null
            ? XKCDColors.SunshineYellow
            : fixed_body.orbitDriver.Renderer.nodeColor;
    } else if (type == MapObject.ObjectType.ApproachIntersect) {
      associated_map_object = reference_frame.target_override.mapObject;
      colour = XKCDColors.Chartreuse;
    } else {
      if (!reference_frame.target_override &&
          (reference_frame.frame_type ==
               ReferenceFrameSelector.FrameType.BODY_CENTRED_NON_ROTATING ||
           reference_frame.frame_type ==
               ReferenceFrameSelector.FrameType.BODY_SURFACE)) {
        // In one-body frames, the apsides are shown with the colour of the
        // body.
        // The nodes are with respect to the equator, rather than with respect
        // to an orbit. We show the nodes in a different (but arbitrary) colour
        // so that they can be distinguished easily.
        associated_map_object = reference_frame.selected_celestial.MapObject;
        colour = XKCDColors.Chartreuse;
      } else {
        // In two-body frames, if apsides are shown, they are shown with the
        // colour of the secondary (or in XKCD chartreuse if the secondary is a
        // vessel).
        // The nodes are with respect to the orbit of the secondary around the
        // primary. We show the nodes with the colour of the primary.
        CelestialBody primary = reference_frame.target_override
            ? reference_frame.selected_celestial
            : reference_frame.selected_celestial.referenceBody;
        associated_map_object = primary.MapObject;
        colour = primary.orbit == null
            ? XKCDColors.SunshineYellow
            : primary.orbitDriver.Renderer.nodeColor;
      }
    }
    colour.a = 1;

    for (int i = 0; i < 64 && !apsis_iterator.IteratorAtEnd();
         ++i, apsis_iterator.IteratorIncrement()) {
      QP apsis = apsis_iterator.IteratorGetDiscreteTrajectoryQP();
      MapNodeProperties node_properties = new MapNodeProperties {
        visible = true,
        object_type = type,
        colour = colour,
        reference_frame = reference_frame,
        world_position = (Vector3d)apsis.q,
        velocity = (Vector3d)apsis.p,
        source = source,
        time = apsis_iterator.IteratorGetDiscreteTrajectoryTime()
      };
      if (type == MapObject.ObjectType.Periapsis &&
          reference_frame.selected_celestial.GetAltitude(
              node_properties.world_position) < 0) {
        node_properties.object_type = MapObject.ObjectType.PatchTransition;
        node_properties.colour = XKCDColors.Orange;
      }

      if (pool_index_ == nodes_.Count) {
        nodes_.Add(MakePoolNode());
      } else if (properties_[nodes_[pool_index_]].object_type !=
                     node_properties.object_type ||
                 properties_[nodes_[pool_index_]].colour !=
                     node_properties.colour) {
        // KSP attaches labels to its map nodes, but never detaches them.
        // If the node changes type, we end up with an arbitrary combination of
        // labels Ap, Pe, AN, DN.
        // Similarly, if the node changes celestial, the colour of the icon
        // label is not updated to match the icon (making it unreadable in some
        // cases). Recreating the node entirely takes a long time (approximately
        // 𝑁 * 70 μs, where 𝑁 is the total number of map nodes in existence),
        // instead we manually get rid of the labels.
        foreach (var component in
                 nodes_[pool_index_].transform.GetComponentsInChildren<
                     TMPro.TextMeshProUGUI>()) {
          if (component.name == "iconLabel(Clone)") {
            UnityEngine.Object.Destroy(component.gameObject);
          }
        }
        // Ensure that KSP thinks the type changed, and reattaches icon
        // labels next time around, otherwise we might end up with no labels.
        // Null nodes do not have a label, so inducing a type change through
        // Null does not result in spurious labels.  Note that the type is
        // updated only if the node is visible.
        properties_[nodes_[pool_index_]].visible = true;
        properties_[nodes_[pool_index_]].object_type =
            MapObject.ObjectType.Null;
        nodes_[pool_index_].NodeUpdate();
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
            var properties = properties_[node];
            if (PlanetariumCamera.fetch.target !=
                properties_[node].associated_map_object) {
              PlanetariumCamera.fetch.SetTarget(
                  properties_[node].associated_map_object);
            }
          }
        };
    new_node.OnUpdateVisible +=
        (KSP.UI.Screens.Mapview.MapNode node,
         KSP.UI.Screens.Mapview.MapNode.IconData icon) => {
          icon.visible = properties_[node].visible;
          icon.color = properties_[node].colour;
        };
    new_node.OnUpdateType +=
        (KSP.UI.Screens.Mapview.MapNode node,
         KSP.UI.Screens.Mapview.MapNode.TypeData type) => {
          var properties = properties_[node];
          type.oType = properties.object_type;
          if (properties.object_type ==
              MapObject.ObjectType.PatchTransition) {
            type.pType =
                KSP.UI.Screens.Mapview.MapNode.PatchTransitionNodeType.Impact;
          } else if (properties.object_type ==
                     MapObject.ObjectType.ApproachIntersect) {
            type.aType = KSP.UI.Screens.Mapview.MapNode.ApproachNodeType
                             .CloseApproachOwn;
          }
        };
    new_node.OnUpdateCaption +=
        (KSP.UI.Screens.Mapview.MapNode node,
         KSP.UI.Screens.Mapview.MapNode.CaptionData caption) => {
          var properties = properties_[node];
          string source;
          switch (properties.source) {
            case NodeSource.FLIGHT_PLAN:
              source = "Planned";
              break;
            case NodeSource.PREDICTION:
              source = "Predicted";
              break;
            default:
              throw Log.Fatal($"Unexpected node source {properties.source}");
          }
          switch (properties.object_type) {
            case MapObject.ObjectType.Periapsis:
            case MapObject.ObjectType.Apoapsis: {
              String apsis_name =
                  properties.object_type == MapObject.ObjectType.Periapsis
                      ? "Periapsis"
                      : "Apoapsis";
              CelestialBody celestial =
                  properties.reference_frame.selected_celestial;
              Vector3d position = properties.world_position;
              double speed = properties.velocity.magnitude;
              caption.Header =
                  $@"{source} {celestial.name} {apsis_name} :\n{
                     celestial.GetAltitude(position):N0} m".ToString(
                      Culture.culture);
              caption.captionLine2 =
                  $"{speed:N0} m/s".ToString(Culture.culture);
              break;
            }
            case MapObject.ObjectType.AscendingNode:
            case MapObject.ObjectType.DescendingNode: {
              string node_name =
                  properties.object_type == MapObject.ObjectType.AscendingNode
                      ? "Ascending Node"
                      : "Descending Node";
              string plane =
                  properties.reference_frame.ReferencePlaneDescription();
              caption.Header = $"{source} {node_name} :\n{plane}";
              caption.captionLine2 =
                  $"{properties.velocity.z:N0} m/s".ToString(Culture.culture);
              break;
            }
            case MapObject.ObjectType.ApproachIntersect: {
              Vessel target_vessel = properties.reference_frame.target_override;
              double separation = (target_vessel.GetWorldPos3D() -
                                   properties.world_position).magnitude;
              double speed = properties.velocity.magnitude;
              caption.Header =
                  $@"{source} Target Approach : {separation:N0} m".ToString(
                      Culture.culture);
              caption.captionLine2 =
                  $"{speed:N0} m/s".ToString(Culture.culture);
              break;
            }
            case MapObject.ObjectType.PatchTransition: {
              CelestialBody celestial =
                  properties.reference_frame.selected_celestial;
              caption.Header = $"{source} {celestial.name} Impact";
              caption.captionLine1 = "";
              caption.captionLine2 = "";
              break;
            }
          }
          if (properties.object_type != MapObject.ObjectType.PatchTransition) {
            caption.captionLine1 =
                "T" + FlightPlanner.FormatTimeSpan(TimeSpan.FromSeconds(
                          Planetarium.GetUniversalTime() - properties.time));
          }
        };
    new_node.OnUpdatePosition +=
        (KSP.UI.Screens.Mapview.MapNode node) =>
        ScaledSpace.LocalToScaledSpace(properties_[node].world_position);
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
    public UnityEngine.Color colour;
    public MapObject associated_map_object;
    public ReferenceFrameSelector reference_frame;
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
