# This scripts takes the output of the HORIZONS service, with the options used
# in query.f, and produces a protocol buffer text format file for a
# |principia.serialization.SolarSystemFile| containing an |initial_state| field.
BEGIN {
  print "initial_state {";
  print "  needs : \"RealSolarSystem\"";
  print "  frame : ICRF_J2000_EQUATOR";
  n = 0;
}
/^Target body name:/ {
  body = $0;
  sub(/ +\{.*/, "", body);
  sub(/ \(.*/, "", body);
  sub(/.*: ([0-9]+ )?/, "", body)
  bodies[n++] = body;
}
/^[0-9]+\.[0-9]+ = .* \(CT\)/ {
  if (jd == "") {
    jd = $1;
  } else if (jd != $1) {
    print "ERROR: inconsistent dates";
    exit 1;
  }
}
/^ *X = / {
  x[body] = $3;
  y[body] = $6;
  z[body] = $9
}
/^ *VX=/ {
  vx[body] = $2;
  vy[body] = $4;
  vz[body] = $6;
}
END {
  print "  # The time of the initial state, as a TDB Julian date.";
  if (jd == "2433282.500000000") {
    print "  # This is 1950-01-01T00:00:00 TDB.";
    print "  epoch : " jd
  } else {
    print "ERROR: unexpected epoch";
    exit 1
  }

  # The year is 2015 and I am writing a Bubble Sort.  I kid you not.
  do {
    swapped = 0;
    for (i = 0; i < n; ++i) {
      if (bodies[i - 1] > bodies[i]) {
        t = bodies[i]
        bodies[i] = bodies[i - 1]
        bodies[i - 1] = t
        swapped = 1
      }
    }
  } while (swapped);

  for (i = 0; i < n; ++i) {
    b = bodies[i];
    print "  body {";
    print "    name : \"" b "\"";
    print "    x    : \"" x[b] " km\"";
    print "    y    : \"" y[b] " km\"";
    print "    z    : \"" z[b] " km\"";
    print "    vx   : \"" vx[b] " km/s\"";
    print "    vy   : \"" vy[b] " km/s\"";
    print "    vz   : \"" vz[b] " km/s\"";
    print "  }"
  }

  print "}"
}
