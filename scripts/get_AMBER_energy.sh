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

# The script parses the output of NAB of AmberTools and greps the energy terms. 
# The input file is taken as a cmdline argument.

fn=$1

if [ -z "$fn" ]; then
    echo >&2 -e "--- !Error in $0: no cmd line arguments given. Aborting \n  Usage: $0 naboutfilename \n Aborting.\n"; exit 1
fi

Total=`cat ${fn} | grep -m 1  "ff:" | awk '{print $3}'`
bad=`cat ${fn} | grep -m 1  "ff:" | awk '{print $4}'`
vdW=`cat ${fn} | grep -m 1  "ff:" | awk '{print $5}'`
elect=`cat ${fn} | grep -m 1  "ff:" | awk '{print $6}'`

fact=`unitconv.py kcal/mol toGMX`
python -c "print 'E_POT '+ str($Total * $fact )"
python -c "print 'E_BONDED '+ str( $bad * $fact )"
python -c "print 'E_VDW '+ str( $vdW * $fact )"
python -c "print 'E_ELEC '+ str( $elect * $fact )"

