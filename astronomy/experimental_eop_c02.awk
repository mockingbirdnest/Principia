BEGIN {
  out = ""; 
  entries = 0
}
$1 != "!" {
  ++entries;
  out = out sprintf("  ExperimentalEOPC02Entry(%10.3f, %11.7f * Second),\n",
                    $1, $4)
}
END {
  printf("constexpr std::array<ExperimentalEOPC02Entry, %d> " \
             "experimental_eop_c02 = {{\n",
         entries);
  printf(out);
  print "}};"
}
