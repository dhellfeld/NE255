#!/bin/sh

# sh macrowriter.sh <nodes> <particles> <nside>

rm run*

nodes=$1
particles=$2
nside=$3

hpi=$((12*$nside*$nside))

for i in `seq 1 $nodes`; do

echo "
# Initialize kernel
/run/initialize

/control/verbose 0
/tracking/verbose 0
/run/verbose 0
/event/verbose 0

/PRISM_SIM/geom/fullMask

/PRISM_SIM/run/DOI true

# Configure HEALPix
/PRISM_SIM/direction/setHPNside $nside
/PRISM_SIM/direction/setHPindexscheme Ring
/PRISM_SIM/direction/SetUpHEALPix

/PRISM_SIM/run/particles $particles
/PRISM_SIM/run/HPstart $((($i - 1) * ($hpi / $nodes) + 1))

/run/beamOn $(($hpi * $particles / $nodes))
" >> run$i.mac;

done

