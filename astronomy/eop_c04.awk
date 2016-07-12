BEGIN {
  out = ""; 
  entries = 0
}
start_of_series || $1 == 1962 {
  start_of_series = 1
}
start_of_series {
  ++entries;
  out = out sprintf("  UTCToUT1MinusUTC(\"%04d-%02d-%02dT00:00:00\"_DateTime, %10.7f * Second),\n", $1, $2, $3, $7)
}
END {
  printf("constexpr std::array<UTCToUT1MinusUTC, %d> eop_c04 = {{\n", entries);
  printf(out);
  print "}};"
}