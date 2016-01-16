BEGIN {
  print "["
}

/^(.*)_body\.hpp/ {
  body = $0
  header = body
  sub(/_body\.hpp$/, ".hpp", header)
  print "{ include: [\"\\\"" body "\\\"\", private, \"\\\"" header "\\\"\", public] },"
}

END {
  print "]"
}
