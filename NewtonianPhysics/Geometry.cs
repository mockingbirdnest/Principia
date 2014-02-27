namespace NewtonianPhysics.Geometry {
  // Formalism for duration, length, speed and acceleration.

  #region Dimensions

  public struct Acceleration {
    public double value;
    public static Acceleration operator -(Acceleration a1, Acceleration a2) {
      return new Acceleration { value = a1.value - a2.value };
    }
    public static Acceleration operator -(Acceleration a) {
      return new Acceleration { value = -a.value };
    }
    public static Acceleration operator +(Acceleration a1, Acceleration a2) {
      return new Acceleration { value = a1.value + a2.value };
    }
  }
  public struct Duration {
    public double value;
    public static Duration operator -(Duration t1, Duration t2) {
      return new Duration { value = t1.value - t2.value };
    }
    public static Duration operator -(Duration t) {
      return new Duration { value = -t.value };
    }
    public static Speed operator *(Acceleration a, Duration t) {
      return new Speed { value = a.value * t.value };
    }
    public static Speed operator *(Duration t, Acceleration a) {
      return new Speed { value = t.value * a.value };
    }
    public static Length operator *(Duration t, Speed v) {
      return new Length { value = v.value * t.value };
    }
    public static Length operator *(Speed v, Duration t) {
      return new Length { value = t.value * v.value };
    }
    public static Speed operator /(Length d, Duration t) {
      return new Speed { value = d.value / t.value };
    }
    public static Acceleration operator /(Speed v, Duration t) {
      return new Acceleration { value = v.value / t.value };
    }
    public static Duration operator +(Duration t1, Duration t2) {
      return new Duration { value = t1.value + t2.value };
    }
    public explicit operator double(Duration t) {
      return t.value;
    }
    public explicit operator Duration(double t) {
      return new Duration { value = t };
    }
  }
  public struct Length {
    public double value;
    public static Length operator -(Length d1, Length d2) {
      return new Length { value = d1.value - d2.value };
    }
    public static Length operator -(Length d) {
      return new Length { value = -d.value };
    }
    public static Length operator +(Length d1, Length d2) {
      return new Length { value = d1.value + d2.value };
    }
    public explicit operator double(Length d) {
      return d.value;
    }
    public explicit operator Length(double d) {
      return new Length { value = d };
    }
  }
  public struct Speed {
    public double value;
    public static Speed operator -(Speed v1, Speed v2) {
      return new Speed { value = v1.value - v2.value };
    }
    public static Speed operator -(Speed v) {
      return new Speed { value = -v.value };
    }
    public static Speed operator +(Speed v1, Speed v2) {
      return new Speed { value = v1.value + v2.value };
    }
    public explicit operator double(Speed v) {
      return v.value;
    }
    public explicit operator Speed(double v) {
      return new Speed { value = v };
    }
  }

  #endregion Dimensions

  // Formalism for Galilean reference frames. Duration does not change between
  // reference frames so there is no need for Duration<T> (this would be needed
  // for special relativity).

  #region Reference Frames

  public struct AbsolutePosition<T> where T : ReferenceFrame<T> {
    public Length x, y, z;
    public static RelativePosition<T> operator -(AbsolutePosition<T> q1,
                                                 AbsolutePosition<T> q2) {
      return new RelativePosition<T> {
        x = q1.x - q2.x,
        y = q1.y - q2.y,
        z = q1.z - q2.z
      };
    }
  }
  public struct AbsoluteTime<T> where T : ReferenceFrame<T> {
    public Duration t;
    public static Duration operator -(AbsoluteTime<T> t1, AbsoluteTime<T> t2) {
      return t1.t - t2.t;
    }
  }
  public struct AbsoluteVelocity<T> where T : ReferenceFrame<T> {
    public Speed x, y, z;
    public static RelativeVelocity<T> operator -(AbsoluteVelocity<T> v1,
                                                 AbsoluteVelocity<T> v2) {
      return new RelativeVelocity<T> {
        x = v1.x - v2.x,
        y = v1.y - v2.y,
        z = v1.z - v2.z
      };
    }
  }
  public struct Acceleration<T> where T : ReferenceFrame<T> {
    public Acceleration x, y, z;
    public static Acceleration<T> operator -(Acceleration<T> a) {
      return new Acceleration<T> { x = -a.x, y = -a.y, z = -a.z };
    }
    public static Acceleration<T> operator -(Acceleration<T> a1,
                                                Acceleration<T> a2) {
      return new Acceleration<T> {
        x = a1.x - a2.x,
        y = a1.y - a2.y,
        z = a1.z - a2.z
      };
    }
    public static Acceleration<T> operator +(Acceleration<T> a1,
                                          Acceleration<T> a2) {
      return new Acceleration<T> {
        x = a1.x + a2.x,
        y = a1.y + a2.y,
        z = a1.z + a2.z
      };
    }
  }
  public struct Event<T> where T : ReferenceFrame<T> {
    public AbsolutePosition<T> q;
    public AbsoluteTime<T> t;
  }
  public struct RelativePosition<T> where T : ReferenceFrame<T> {
    public Length x, y, z;
    public static RelativePosition<T> operator -(RelativePosition<T> q) {
      return new RelativePosition<T> { x = -q.x, y = -q.y, z = -q.z };
    }
    public static RelativePosition<T> operator -(RelativePosition<T> q1,
                                      RelativePosition<T> q2) {
      return new RelativePosition<T> {
        x = q1.x - q2.x,
        y = q1.y - q2.y,
        z = q1.z - q2.z
      };
    }
    public static RelativePosition<T> operator +(RelativePosition<T> q1,
                                          RelativePosition<T> q2) {
      return new RelativePosition<T> {
        x = q1.x + q2.x,
        y = q1.y + q2.y,
        z = q1.z + q2.z
      };
    }
    public static AbsolutePosition<T> operator +(AbsolutePosition<T> q1,
                                                 RelativePosition<T> q2) {
      return new AbsolutePosition<T> {
        x = q1.x + q2.x,
        y = q1.y + q2.y,
        z = q1.z + q2.z
      };
    }
    public static AbsolutePosition<T> operator +(RelativePosition<T> q1,
                                                 AbsolutePosition<T> q2) {
      return q2 + q1;
    }
  }
  public struct RelativeVelocity<T> where T : ReferenceFrame<T> {
    public Speed x, y, z;
    public static RelativeVelocity<T> operator -(RelativeVelocity<T> v) {
      return new RelativeVelocity<T> { x = -v.x, y = -v.y, z = -v.z };
    }
    public static RelativeVelocity<T> operator -(RelativeVelocity<T> v1,
                                                RelativeVelocity<T> v2) {
      return new RelativeVelocity<T> {
        x = v1.x - v2.x,
        y = v1.y - v2.y,
        z = v1.z - v2.z
      };
    }
    public static RelativeVelocity<T> operator +(RelativeVelocity<T> v1,
                                          RelativeVelocity<T> v2) {
      return new RelativeVelocity<T> {
        x = v1.x + v2.x,
        y = v1.y + v2.y,
        z = v1.z + v2.z
      };
    }
    public static AbsoluteVelocity<T> operator +(AbsoluteVelocity<T> v1,
                                                 RelativeVelocity<T> q2) {
      return new AbsoluteVelocity<T> {
        x = v1.x + q2.x,
        y = v1.y + q2.y,
        z = v1.z + q2.z
      };
    }
    public static AbsoluteVelocity<T> operator +(RelativeVelocity<T> q1,
                                                 AbsoluteVelocity<T> q2) {
      return q2 + q1;
    }
  }
  public abstract class ReferenceFrame<T> where T : ReferenceFrame<T> {
    public Event<T> origin = new Event<T> {
      q = new AbsolutePosition<T> {
        x = (Length)0,
        y = (Length)0,
        z = (Length)0
      },
      t=new AbsoluteTime<T> {t
        =(Duration)0}
    };
  }

  #endregion Reference Frames

  #region Frame Types

  public abstract class InertialReferenceFrame<T> : ReferenceFrame<T>
  where T : InertialReferenceFrame<T> { };
  public abstract class UniformlyRotatingReferenceFrame<T> :
  ReferenceFrame<T>
   where T : UniformlyRotatingReferenceFrame<T> {
    public AngularVelocity<T> ω;
    public explicit  operator Event<InertialFrame> (Event<T> e){
      return new Event<InertialFrame> {
        q= ω * (e.t-origin.t) * AbsolutePosition<>
      }
    }
public class InertialFrame:InertialReferenceFrame<InertialFrame> {
    }
  };

  #endregion Frame Types

#region Rotations

    public struct AngularVelocity<T> where T:ReferenceFrame<T>  {
          public SpatialRotation<T> operator*(AngularVelocity<T> ω, Duration t){
      return new SpatialRotation<T>{};
    }
  }
public struct SpatialRotation<T> where T:ReferenceFrame<T> {
   private Versor rotation;
    public AbsolutePosition<T> operator*(SpatialRotation<T> r, AbsolutePosition<T>){
      return new AbsolutePosition<T>{};
    }
  }

#endregion Rotations

#region Pure Geometry

  public struct PseudoVector{
  }
public struct Vector{
  }
public struct Versor{
  }

#endregion Pure Geometry
}