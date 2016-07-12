BEGIN {
  out = ""; 
  entries = 0
}
$1 != "!" {
  ++entries;
  out = out sprintf("  MJDToUT1MinusTAI(%10.3f, %11.7f * Second),\n",
                    $1, $4)
}
END {
  printf("constexpr std::array<MJDToUT1MinusTAI, %d> " \
             "const experimental_eop_c02 = {{\n",
         entries);
  printf(out);
  print "}};"
}
