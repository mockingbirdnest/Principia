curl 'https://hpiers.obspm.fr/iers/series/longterm/eopc02.1830-now' |
  sed 's/-/ -/g' |
  awk -f experimental_eop_c02.awk  > experimental_eop_c02.generated.h

curl 'https://hpiers.obspm.fr/iers/eop/eopc04/eopc04.62-now' |
  awk -f eop_c04.awk > eop_c04.generated.h
