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


###------------------------------------------------------------------------------
### Default parameters

pp=./confs
FFCONV_PATH=`which ffconv.py | sed "s:/ffconv.py::g"`

d=FFLDSERVER

BriefUsage="--------------------------------------------------------------------------
 BRIEF USAGE (for $0)

    $0 --mol paracetamol --confpath ./confs --ffconvpath ~/ffconv

 OPTIONS
    --mol            - name of the molecule (defines the names of input files)
    --confpath       - [default: ${confpath}] path to where PDB file(s) with coordinates of molecule located.
                       if many present and lskipcheck != yes: then energies of configuration from all pdb files 
                       are used to check the correctness of conversion.
    --ffconvpath     - [default: \$FFCONVPATH] path to ffconv distribution. If FFCONV not defined then the path set 
                       to location of ffconv.py executable.
    --lskipcheck     - skip energy evaluation
    -h               - print help

"
HelpMessage="
HELP MESSAGE FOR $0

The script performs a simple enegy evaluation with CHARMM.

${BriefUsage}

"


if [ -z "$*" ]; then
    echo >&2 -e "$BriefUsage"
    echo >&2 -e "! Error in $0: no cmd line arguments given. Aborting.\n"; exit 1;
fi

###------------------------------------------------------------------------------
### Parse command line arguments
echo -e "Running: $0 $* \n"
while [ $# -gt 0 ] ; do
  case $1
  in
    --confpath)
      pp=$2; shift 2
    ;;
    --mol)
      mol=$2; shift 2 
    ;;
    --ffconvpath)
      FFCONV_PATH=$2; shift 2 
    ;;
    --lskipcheck)
      lskipcheck=$2; shift 2 
    ;;
    -h)
      shift 1;
      echo -e "$HelpMessage"; exit 0
    ;;
    *)
      echo >&2 -e "$BriefUsage"
      echo >&2 -e "!Error in "$0": the command line argument [$1] not known. Aborting.\n" ; exit 1;
    ;;
  esac
done

###------------------------------------------------------------------------------
### Checks

if [ -z "$mol" ]; then
      echo >&2 -e "--- !Error in "$0": molecule name NOT set. Use --mol key. Aborting.\n" ; exit 1;
fi


###------------------------------------------------------------------------------
### Set path
templdir=$FFCONV_PATH/templates/$d
scriptsdir=$FFCONV_PATH/scripts

if [ "$lskipcheck" != 'yes' ]; then

if [ -z "$SCHRODINGER" ]; then
      echo >&2 -e "--- !Error in "$0": SCHRODINGER environmental variable NOT set. Aborting.\n" ; exit 1;
fi


###------------------------------------------------------------------------------
### Do energy evaluation for each configuration 
mkdir $d
for i in `ls ${pp}/*.pdb`; do
    conf=`python -c "t=\"$i\".split('/'); t=t[-1]; t=t.split('.'); t=t[0]; print t"`

    mkdir $d/$conf 
    p=`pwd`; cd $d/$conf

        # Energy evaluation
        cmd="$SCHRODINGER/utilities/ffld_server -ipdb $i -print_energy_components -version 2005 > energy.out"
        eval $cmd ; if [ $? -ne 0 ]; then { echo >&2 -e "--- !Error: [$cmd] failed. Aborting.\n" ; exit 1; } fi

        # Grep energy data to E.dat file 
        $scriptsdir/get_FFLDSERVER_energy.sh energy.out > E.dat

    cd $p
done


fi # lskipcheck


### Files to store to DB


if [ ! -f "${mol}.pdb" ]; then
    echo >&2 -e "--- !Note in ["$0"]: no ${mol}.pdb present. Copy the last entry from the confpath. Continue.\n" ;
    cp ${pp}/${conf}.pdb ${mol}.pdb
fi
    

echo "${mol}.pdb
${mol}.ffld
" >> MOLFILES


