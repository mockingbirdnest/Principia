BEGIN {
  print "principia_initial_state {";
  first_pass = true
}
/^Target body name:/ {
  if (first_pass) {
    first_pass = false
  } else {
    print "  body {";
    print "    name = " body;
    print "    x    = " x " km";
    print "    y    = " y " km";
    print "    z    = " z " km";
    print "    vx   = " vx " km/s";
    print "    vy   = " vy " km/s";
    print "    vz   = " vz "km/s";
    print "  }"
  }
  body = $0;
  sub(/ \(.*/, "", body);
  sub(/.*: /, "", body)
}
/^ *X = / {
  x = $3;
  y = $6;
  z = $9
}
/^ *VX=/ {
  vx = $2;
  vy = $4;
  vz = $6
}
END {
  print "}"
}
