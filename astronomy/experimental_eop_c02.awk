# Processes the Experimental EOP C02 time series
# https://hpiers.obspm.fr/iers/series/longterm/eopc02.1830-now, see also
# https://hpiers.obspm.fr/iers/series/longterm/eopc02.txt.
# Field 1 is the MJD, field 4 is UT1 - TAI in seconds.
BEGIN {
  out = ""; 
}
$1 != "#" {
  out = out sprintf("  ExperimentalEOPC02Entry(%10.3f, %11.7f * Second),\n",
                    $1, $4)
}
END {
  printf("constexpr std::array experimental_eop_c02{\n");
  printf(out);
  print "};"
}
