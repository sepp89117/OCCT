# Synthetic periodic seam fillet fixture. No user model data is used.
# Full cylinder R20 x H10, quarter-cylinder R20 x H10 above it.
# Rotate the lower cylinder so its technical seam meets one radial edge.
pload MODELING
pcylinder base 20 10
trotate base 0 0 0 0 0 1 90
pcylinder boss 20 10 90
ttranslate boss 0 0 10
bfuse body base boss
checkshape body
bopcheck body
save body PeriodicSeamSynthetic.brep
exit
