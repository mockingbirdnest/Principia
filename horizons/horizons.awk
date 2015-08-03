BEGIN {
  print "principia_initial_state {"
}
/^Target body name:/ {
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
  print "}"
}
