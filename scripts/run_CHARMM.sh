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
CHARMM="charmm"

pp=./confs
FFCONV_PATH=`which ffconv.py | sed "s:/ffconv.py::g"`

prog=CHARMM

d=$prog

BriefUsage="--------------------------------------------------------------------------
BRIEF USAGE (for $0)

OPTIONS

    --mol            - name of the molecule (defines the names of input files)
    --confpath       - [default: ${confpath}] path to where PDB file(s) with coordinates of molecule located.
                       if many present and lskipcheck != yes: then energies of configuration from all pdb files 
                       are used to check the correctness of conversion.
    --ffconvpath     - [default: \$FFCONVPATH] path to ffconv distribution. If FFCONV not defined then the path set 
                       to location of ffconv.py executable.
    --charmmcmd      - [default: $CHARMM] executable for namd
    --charmm-ff-prm  - charmm PRM force field file (only used when to- or fromprog is CHARMM)
    --charmm-ff-rtf  - charmm RTF force field file (only used when to- or fromprog is CHARMM)
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
    --charmmcmd)
      CHARMM=$2; shift 2 
    ;;
    --charmm-ff-prm)
      charmmffprm=$2; shift 2 
    ;;
    --charmm-ff-rtf)
      charmmffrtf=$2; shift 2 
    ;;
#    --fstr)
#      fstr=$2; shift 2 
#    ;;
    --lskipcheck)
      lskipcheck=$2; shift 2 
    ;;
    -h)
      shift 1;
      echo -e "$HelpMessage"; exit 0
    ;;
    *)
      echo >&2 -e "$BriefUsage"
      echo >&2 -e "--- !Error in [ "$0" ]: the command line argument [$1] not known. Aborting.\n" ; exit 1;
    ;;
  esac
done




### Vars
templdir=$FFCONV_PATH/templates/$prog
scriptsdir=$FFCONV_PATH/scripts

if [ "$lskipcheck" != 'yes' ]; then

mkdir $d

for i in `ls $pp/*.pdb`; do
    conf=`python -c "t=\"$i\".split('/'); t=t[-1]; t=t.split('.'); t=t[0]; print t"`

    mkdir $d/$conf 
    p=`pwd`; cd $d/$conf

        fpdb=$i
        fstr=../../${mol}.str
        cmol=`python -c "t='$mol'; print t[0:4].upper()"`
        FillVariablesToTemplate.sh FONAME=charmm.in  TEMPLATE=${templdir}/template_charmm.in TMPPDB=${fpdb} TMPSTR=${fstr} TMPPRM=${charmmffprm} TMPRTF=${charmmffrtf}  TMPMOL=${cmol}

        cmd="$CHARMM < charmm.in > charmm.out"
        eval $cmd ; if [ $? -ne 0 ]; then { echo >&2 -e "--- !Error: [$cmd] failed. Aborting.\n" ; exit 1; } fi

        lfailed=`tail -n 20 charmm.out | grep "ABNORMAL TERMINATION" | wc -l`
        if [ "$lfailed" -ne 0 ]; then
             echo >&2 -e "--- !Error: [$cmd] failed. Charmm has terminated abnormally. Aborting.\n" ; exit 1; 
        fi
        
        $scriptsdir/get_CHARMM_energy.py charmm.out > E.dat

    cd $p
done


fi # lskipcheck



### Files to store to DB

if [ ! -f "${mol}.pdb" ]; then
    echo >&2 -e "--- !Note in ["$0"]: no ${mol}.pdb present. Copy the last entry from the confpath. Continue.\n" ;
    cp ${pp}/${conf}.pdb ${mol}.pdb
fi
    

echo "${mol}.pdb
${mol}.str
" >> MOLFILES

echo "${charmmffprm}
${charmmffrtf}
" >> FFFILES



