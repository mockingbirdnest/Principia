BEGIN {
  out = ""; 
  entries = 0
}
start_of_series || $1 == 1962 {
  start_of_series = 1
}
start_of_series {
  ++entries;
  out = out sprintf("  EOPC04Entry(%d'%02d'%02d, %10.7f * Second),\n",
                    $1, $2, $3, $7)
}
END {
  printf("constexpr std::array<EOPC04Entry, %d> eop_c04 = {{\n", entries);
  printf(out);
  print "}};"
}
