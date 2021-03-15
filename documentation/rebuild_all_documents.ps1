# These files are unmaintained; they may not compile.
$ancient_files = @('Hamiltonian Mechanics', 'ODEs and Runge-Kutta integrators')

# The UTC second nearest J2000 (which would be 946727935.816 using TT(TAI)).
$env:SOURCE_DATE_EPOCH=946727936

foreach ($f in ls *.tex) {
  for ($i = 0; $i -lt 3; ++$i) {
    if ($ancient_files.contains($f.basename)) {
      continue
    }
    if (test-path "$($f.basename).bcf") {
     biber $f.basename
     if (-not $?) { return $? }
    }
    xelatex $f.basename
    if (-not $?) { return $? }
  }
}
