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

# This script parses the ffld_server output file and greps the energy terms.
# Input file is taken as a cmdline argument.

fn=$1

if [ -z "$fn" ]; then
    echo >&2 -e "--- !Error in $0: no cmd line arguments given. Aborting \n  Usage: $0 ffldoutfilename \n Aborting.\n"; exit 1
fi


#Total Energy             -5.949611 
#Bond Stretch              1.821006 
#Angle Bending             2.669799 
#Torsion Angle             2.062473 
#Improper Torsion          0.016145 
#1,4 Lennard Jones         7.762752 
#1,4 Electrostatic       -39.489182 
#Lennard Jones             1.811094 
#Electrostatic            17.396301 
#Torsion Restraints        0.000000 

Total=`cat ${fn} | grep   "Total Energy" | awk '{print $3}'`
Bond=`cat ${fn} | grep   "Bond Stretch" | awk '{print $3}'`
Angle=`cat ${fn} | grep   "Angle Bending" | awk '{print $3}'`
Torsion=`cat ${fn} | grep   "Torsion Angle" | awk '{print $3}'`
Improper=`cat ${fn} | grep  "Improper Torsion" | awk '{print $3}'`
LJ14=`cat ${fn} | grep  "1,4 Lennard Jones" | awk '{print $4}'`
Elec14=`cat ${fn} | grep   "1,4 Electrostatic" | awk '{print $3}'`
LJ=`cat ${fn} | grep   "^Lennard Jones " | awk '{print $3}'`
Elec=`cat ${fn} | grep   "^Electrostatic " | awk '{print $2}'`
Rest=`cat ${fn} | grep   "Torsion Restraints" | awk '{print $3}'`

fact=`unitconv.py kcal/mol toGMX`
python -c "print 'E_POT '+ str($Total * $fact )"
python -c "print 'E_BONDS '+ str($Bond * $fact )"
python -c "print 'E_ANGLES '+ str($Angle * $fact )"
python -c "print 'E_TORSIONS '+ str($Torsion * $fact )"
python -c "print 'E_IMPROPERS '+ str($Improper * $fact )"
python -c "print 'E_DIHEDRALS '+ str( ( $Torsion + $Improper ) * $fact )"
python -c "print 'E_VDW '+ str( ( $LJ + $LJ14 ) * $fact )"
python -c "print 'E_ELEC '+ str( ( $Elec + $Elec14) * $fact )"
python -c "print 'E_BONDED '+ str( ( $Bond + $Angle + $Torsion + $Improper ) * $fact )"




