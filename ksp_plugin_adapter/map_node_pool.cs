using System;
using System.Collections.Generic;
using KSP.Localization;

namespace principia {
namespace ksp_plugin_adapter {

internal class MapNodePool {
  // We render at most 64 markers of one type and one provenance (e.g., at
  // most 64 perilunes for the prediction of the active vessel).  This is
  // more than is readable, and keeps the size of the map node pool under
  // control.
  public const int MaxNodesPerProvenance = 64;

  public enum NodeSource {
    Prediction,
    FlightPlan,
  }

  public struct Provenance {
    public string vessel_guid;
    public NodeSource source;
    public MapObject.ObjectType type;

    public Provenance(string vessel_guid,
                      NodeSource source,
                      MapObject.ObjectType type) {
      this.vessel_guid = vessel_guid;
      this.source = source;
      this.type = type;
    }

    public override int GetHashCode() =>
    (vessel_guid, source, type).GetHashCode();
    public override bool Equals(object other) =>
        other is Provenance p && Equals(p);
    public bool Equals(Provenance other) => vessel_guid == other.vessel_guid &&
                                            source == other.source &&
                                            type == other.type;
  }

  public class SingleProvenancePool {
    public int nodes_used;
    public List<KSP.UI.Screens.Mapview.MapNode> nodes =
        new List<KSP.UI.Screens.Mapview.MapNode>(MaxNodesPerProvenance);
  }

  public MapNodePool(Func<bool> show_only_pinned) {
    nodes_ = new Dictionary<Provenance, SingleProvenancePool>();
    properties_ =
        new Dictionary<KSP.UI.Screens.Mapview.MapNode, MapNodeProperties>();
    show_only_pinned_ = show_only_pinned;
  }

  public void Clear() {
    foreach (var provenance_pool in nodes_) {
      foreach (var node in provenance_pool.Value.nodes) {
        node.Terminate();
      }
    }
    nodes_ = new Dictionary<Provenance, SingleProvenancePool>();
    properties_ =
        new Dictionary<KSP.UI.Screens.Mapview.MapNode, MapNodeProperties>();
  }

  public void Update() {
    List<Provenance> unused_provenances = new List<Provenance>();
    foreach (var provenance_pool in nodes_) {
      var provenance = provenance_pool.Key;
      var pool = provenance_pool.Value;
      if (pool.nodes_used == 0) {
        unused_provenances.Add(provenance);
        continue;
      }
      for (int i = pool.nodes_used; i < pool.nodes.Count; ++i) {
        if (properties_[pool.nodes[i]].visible) {
          properties_[pool.nodes[i]].visible = false;
          pool.nodes[i].NodeUpdate();
        }
      }
      for (int i = 0; i < pool.nodes_used; ++i) {
        pool.nodes[i].NodeUpdate();
      }
      pool.nodes_used = 0;
    }
    foreach (var provenance in unused_provenances) {
      foreach (var node in nodes_[provenance].nodes) {
        node.Terminate();
        properties_.Remove(node);
      }
      nodes_.Remove(provenance);
    }
  }

  public void RenderMarkers(DisposableIterator apsis_iterator,
                            Provenance provenance,
                            ReferenceFrameSelector reference_frame) {
    if (!nodes_.ContainsKey(provenance)) {
      nodes_.Add(provenance, new SingleProvenancePool());
    }
    var pool = nodes_[provenance];
    MapObject associated_map_object;
    UnityEngine.Color colour;
    switch (provenance.type) {
      case MapObject.ObjectType.Apoapsis:
      case MapObject.ObjectType.Periapsis:
        CelestialBody fixed_body = reference_frame.Centre();
        associated_map_object = fixed_body.MapObject;
        colour = fixed_body.orbit == null
                     ? XKCDColors.SunshineYellow
                     : fixed_body.orbitDriver.Renderer.nodeColor;
        break;
      case MapObject.ObjectType.ApproachIntersect:
        associated_map_object = reference_frame.target.mapObject;
        colour = XKCDColors.Chartreuse;
        break;
      case MapObject.ObjectType.AscendingNode:
      case MapObject.ObjectType.DescendingNode:
        if (reference_frame.Centre() == null) {
          // In two-body frames, if apsides are shown, they are shown with the
          // colour of the secondary (or in XKCD chartreuse if the secondary is
          // a vessel).
          // The nodes are with respect to the orbit of the secondary around the
          // primary. We show the nodes with the colour of the primary.
          CelestialBody primary = reference_frame.OrientingBody();
          associated_map_object = primary.MapObject;
          colour = primary.orbit == null
                        ? XKCDColors.SunshineYellow
                        : primary.orbitDriver.Renderer.nodeColor;
        } else {
          // In one-body frames, the apsides are shown with the colour of the
          // body.
          // The nodes are with respect to the equator, rather than with respect
          // to an orbit. We show the nodes in a different (but arbitrary)
          // colour so that they can be distinguished easily.
          associated_map_object = reference_frame.Centre().MapObject;
          colour = XKCDColors.Chartreuse;
        }
        break;
      default:
        throw Log.Fatal($"Unexpected type {provenance.type}");
    }
    colour.a = 1;

    for (int i = 0;
         i < MaxNodesPerProvenance && !apsis_iterator.IteratorAtEnd();
         ++i, apsis_iterator.IteratorIncrement()) {
      QP apsis = apsis_iterator.IteratorGetDiscreteTrajectoryQP();
      MapNodeProperties node_properties = new MapNodeProperties {
          visible = true,
          object_type = provenance.type,
          colour = colour,
          reference_frame = reference_frame,
          world_position = (Vector3d)apsis.q,
          velocity = (Vector3d)apsis.p,
          source = provenance.source,
          time = apsis_iterator.IteratorGetDiscreteTrajectoryTime(),
          associated_map_object = associated_map_object,
      };
      if (provenance.type == MapObject.ObjectType.Periapsis &&
          reference_frame.Centre().GetAltitude(
              node_properties.world_position) < 0) {
        node_properties.object_type = MapObject.ObjectType.PatchTransition;
        node_properties.colour = XKCDColors.Orange;
      }

      if (pool.nodes_used == pool.nodes.Count) {
        pool.nodes.Add(MakePoolNode());
      }
      properties_[pool.nodes[pool.nodes_used++]] = node_properties;
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
        (KSP.UI.Screens.Mapview.MapNode node, Mouse.Buttons buttons) => {
          if (buttons == Mouse.Buttons.Left) {
            MapNodeProperties properties = properties_[node];
            if (PlanetariumCamera.fetch.target !=
                properties.associated_map_object) {
              PlanetariumCamera.fetch.SetTarget(
                  properties.associated_map_object);
            }
          }
        };
    new_node.OnUpdateVisible += (KSP.UI.Screens.Mapview.MapNode node,
                                 KSP.UI.Screens.Mapview.MapNode.IconData
                                     icon) => {
      icon.visible = properties_[node].visible &&
                     (!show_only_pinned_() || node.Pinned);
      icon.color = properties_[node].colour;
    };
    new_node.OnUpdateType += (KSP.UI.Screens.Mapview.MapNode node,
                              KSP.UI.Screens.Mapview.MapNode.TypeData type) => {
      MapNodeProperties properties = properties_[node];
      type.oType = properties.object_type;
      switch (properties.object_type) {
        case MapObject.ObjectType.PatchTransition:
          type.pType = KSP.UI.Screens.Mapview.MapNode.PatchTransitionNodeType.
              Impact;
          break;
        case MapObject.ObjectType.ApproachIntersect:
          type.aType = KSP.UI.Screens.Mapview.MapNode.ApproachNodeType.
              CloseApproachOwn;
          break;
      }
    };
    new_node.OnUpdateCaption += (KSP.UI.Screens.Mapview.MapNode node,
                                 KSP.UI.Screens.Mapview.MapNode.CaptionData
                                     caption) => {
      var properties = properties_[node];
      string source;
      switch (properties.source) {
        case NodeSource.FlightPlan:
          source = L10N.CacheFormat("#Principia_MapNode_Planned");
          break;
        case NodeSource.Prediction:
          source = L10N.CacheFormat("#Principia_MapNode_Predicted");
          break;
        default:
          throw Log.Fatal($"Unexpected node source {properties.source}");
      }
      switch (properties.object_type) {
        case MapObject.ObjectType.Periapsis:
        case MapObject.ObjectType.Apoapsis: {
          CelestialBody celestial = properties.reference_frame.Centre();
          Vector3d position = properties.world_position;
          double speed = properties.velocity.magnitude;
          caption.Header = L10N.CelestialString(
              properties.object_type == MapObject.ObjectType.Periapsis
                  ? "#Principia_MapNode_PeriapsisHeader"
                  : "#Principia_MapNode_ApoapsisHeader",
              new[]{celestial},
              source,
              celestial.GetAltitude(position).FormatN(0));
          caption.captionLine2 =
              L10N.CacheFormat("#Principia_MapNode_ApsisCaptionLine2",
                               speed.FormatN(0));
          break;
        }
        case MapObject.ObjectType.AscendingNode:
        case MapObject.ObjectType.DescendingNode: {
          string node_name =
              properties.object_type == MapObject.ObjectType.AscendingNode
                  ? L10N.CacheFormat("#Principia_MapNode_AscendingNode")
                  : L10N.CacheFormat("#Principia_MapNode_DescendingNode");
          string plane = properties.reference_frame.ReferencePlaneDescription();
          caption.Header = L10N.CacheFormat("#Principia_MapNode_NodeHeader",
                                            source,
                                            node_name,
                                            plane);
          caption.captionLine2 = L10N.CacheFormat(
              "#Principia_MapNode_NodeCaptionLine2",
              properties.velocity.z.FormatN(0));
          break;
        }
        case MapObject.ObjectType.ApproachIntersect: {
          double separation = 
              (properties.reference_frame.target.GetWorldPos3D() -
               properties.world_position).magnitude;
          double speed = properties.velocity.magnitude;
          caption.Header = L10N.CacheFormat("#Principia_MapNode_ApproachHeader",
                                            source,
                                            separation.FormatN(0));
          caption.captionLine2 = L10N.CacheFormat(
              "#Principia_MapNode_ApproachCaptionLine2",
              speed.FormatN(0));
          break;
        }
        case MapObject.ObjectType.PatchTransition: {
          CelestialBody celestial = properties.reference_frame.Centre();
          caption.Header = L10N.CacheFormat("#Principia_MapNode_ImpactHeader",
                                            source,
                                            celestial.Name());
          caption.captionLine1 = "";
          caption.captionLine2 = "";
          break;
        }
      }
      if (properties.object_type != MapObject.ObjectType.PatchTransition) {
        // The font used by map nodes does not support the minus sign, so we
        // fall back to the hyphen-minus.
        caption.captionLine1 =
            "T" + new PrincipiaTimeSpan(
                    Planetarium.GetUniversalTime() - properties.time).
                Format(with_leading_zeroes: false, with_seconds: true)
                .Replace(Culture.culture.NumberFormat.NegativeSign,
                         "-");
      }
    };
    new_node.OnUpdatePosition += (KSP.UI.Screens.Mapview.MapNode node) =>
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

  private Dictionary<Provenance, SingleProvenancePool> nodes_;
  private Dictionary<KSP.UI.Screens.Mapview.MapNode, MapNodeProperties>
      properties_;
  private Func<bool> show_only_pinned_;
}

}  // namespace ksp_plugin_adapter
}  // namespace principia
