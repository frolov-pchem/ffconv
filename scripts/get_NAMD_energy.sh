#!/bin/bash
# Andrey I. Frolov, Jan 2014, ISC RAS, Ivanovo, Russia
#
#   Copyright 2014 Andrey I. Frolov
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.

# The script parses the NAMD output file and grep the energy terms
# The input file is taken as a cmdline argument.

fn=$1

if [ -z "$fn" ]; then
    echo >&2 -e "--- !Error in $0: no cmd line arguments given. Aborting \n  Usage: $0 namdoutput \n Aborting.\n"; exit 1
fi


#ETITLE:      TS           BOND          ANGLE          DIHED          IMPRP               ELECT            VDW       BOUNDARY           MISC        KINETIC               TOTAL           TEMP      POTENTIAL         TOTAL3        TEMPAVG
#ENERGY:       0         0.0000         0.0000         0.0000         0.0000            -80.5907        16.7748         0.0000         0.0000         0.0000            -63.8160         0.0000       -63.8160       -63.8160         0.0000

#ETITLE=`cat ${fn} | grep -m 1  "ETITLE:" | sed "s/ETITLE\://g"`
#ENERGY=`cat ${fn} | grep -m 1  "ETITLE:" | sed "s/ENERGY\://g"`

potential=`cat ${fn} | grep -m 1  "^ENERGY:" | awk '{print $14}'`
vdw=`cat ${fn} | grep -m 1  "^ENERGY:" | awk '{print $8}'`
elect=`cat ${fn} | grep -m 1  "^ENERGY:" | awk '{print $7}'`
bond=`cat ${fn} | grep -m 1  "^ENERGY:" | awk '{print $3}'`
angle=`cat ${fn} | grep -m 1  "^ENERGY:" | awk '{print $4}'`
torsion=`cat ${fn} | grep -m 1  "^ENERGY:" | awk '{print $5}'`
improper=`cat ${fn} | grep -m 1  "^ENERGY:" | awk '{print $6}'`

fact=`unitconv.py kcal/mol toGMX`
python -c "print 'E_POT '+ str($potential * $fact )"
python -c "print 'E_VDW '+ str( $vdw * $fact )"
python -c "print 'E_ELEC '+ str( $elect * $fact )"
python -c "print 'E_BOND '+ str( $bond * $fact )"
python -c "print 'E_ANGLE '+ str( $angle * $fact )"
python -c "print 'E_TORSION '+ str( $torsion * $fact )"
python -c "print 'E_IMPROPER '+ str( $improper * $fact )"
python -c "print 'E_BONDED '+ str( ( $bond + $angle + $torsion + $improper ) * $fact )"

