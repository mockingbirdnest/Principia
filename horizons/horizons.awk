function print_body() {
  print "  body {"
  print "    name = " body;
  print "    x = " x;
  print "    y = " y;
  print "    z = " z;
  print "    vx = " vx;
  print "    vy = " vy;
  print "    vz = " vz;
  print "  }"
}

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
  print_body()
}
END {
  print_body();
  print "}"
}
