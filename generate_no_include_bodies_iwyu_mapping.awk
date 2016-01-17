BEGIN {
  print "["
}

/^(.*)_body\.hpp/ {
  body = $0
  header = body
  sub(/_body\.hpp$/, ".hpp", header)
  if (header != "massless.hpp" &&
      header != "massive.hpp" &&
      header != "oblate.hpp" &&
      header != "rotating.hpp") {
    print "{ include: [\"\\\"" body "\\\"\", private, \"\\\"" header "\\\"\", public] },"
  }
}

END {
  print "]"
}
