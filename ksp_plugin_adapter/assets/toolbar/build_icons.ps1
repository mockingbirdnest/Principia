# To build the icon bitmaps, you need the following installed:
# - [Computer Modern Unicode Fonts](http://cm-unicode.sourceforge.net/download.html)
#   - CMU Serif
#   - CMU Serif Bold
# - [Inkscape](http://inkscape.org)

# Write actual path to Inkscape here:
Set-Alias -name Inkscape `
  -value "C:\Program Files\Inkscape\inkscape.com"

$svgFile = "icons.svg"
$iconSize = 24
$iconNames = @(
  "UI",
  "Settings",
  "PlottingFrame",
  "FlightPlan",
  "TargetCelestial",
  "TargetVessel",
  "PatchedConics",
  "SunLensFlare"
)

foreach($iconIndex in 0..($iconNames.length - 1)) {
  
  $iconName = $iconNames[$iconIndex]
  $x0 = $iconSize * $iconIndex
  $x1 = $x0 + $iconSize
  
  foreach($row in 0..1) {
    
    $state = If ($row) {"on"} Else {"off"}
    $y0 = $iconSize * $row
    $y1 = $y0 + $iconSize
    
    # Making png files
    Inkscape $svgFile `
      --export-png="${iconName}-${state}.png" `
      --export-area=${x0}:${y0}:${x1}:${y1}
  }
}
