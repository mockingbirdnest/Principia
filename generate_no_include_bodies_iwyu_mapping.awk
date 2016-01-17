BEGIN {
  print "["
}

/^(.*)_body\.hpp/ {
  body = $0
  header = body
  sub(/_body\.hpp$/, ".hpp", header)
  if (header != "physics/massless.hpp" &&
      header != "physics/massive.hpp" &&
      header != "physics/oblate.hpp" &&
      header != "physics/rotating.hpp") {
    print "{ include: [\"\\\"" body "\\\"\", private, \"\\\"" header "\\\"\", public] },"
  }
}

END {
  print "]"
}
