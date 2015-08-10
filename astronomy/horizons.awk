# This scripts takes the output of the HORIZONS service, with the options used
# in query.f, and converts produces a KSP config file with nodes named |body|,
# containing values |name|, |x|, |y|, |z|, |vx|, |vy|, |vz|.  |name| is the
# name of the relevant |Celestial| in KSP, the other values are double-precision
# floating point numbers followed by an appropriate unit.
BEGIN {
  print "principia_initial_state {"
}
/^Target body name:/ {
  body = $0;
  sub(/ \(.*/, "", body);
  sub(/.*: ([0-9]+ )?/, "", body)
}
/^[0-9]+\.[0-9]+ = .* \(CT\)/ {
  if (jd == "") {
    jd = $1;
  } else if (jd != $1) {
    print "ERROR: inconsistent dates";
    exit 1
  }
}
/^ *X = / {
  x = $3;
  y = $6;
  z = $9
}
/^ *VX=/ {
  vx = $2;
  vy = $4;
  vz = $6;
  print "  body {";
  print "    name = " body;
  print "    x    = " x " km";
  print "    y    = " y " km";
  print "    z    = " z " km";
  print "    vx   = " vx " km/s";
  print "    vy   = " vy " km/s";
  print "    vz   = " vz " km/s";
  print "  }"
}
END {
  print "  // The time of the initial state, in seconds from the beginning of";
  print "  // the game.  The state in this file is from";
  print "  // JD " jd " TDB.";
  if (jd == "2433282.500000000") {
    print "  // This is 1950-01-01T00:00:00 TDB, we set the epoch one non-leap";
    print "  // year in the past for RSS, so that RSS begins on";
    print "  // 1951-01-01T00:00:00 TDB.";
    print "  epoch = -31536000"
  } else {
    print "ERROR: unexpected epoch";
    exit 1
  }
  print "}"
}
