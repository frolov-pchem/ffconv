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


# TODO : list of keys ?


###------------------------------------------------------------------------------
### Default parameters
NAB="nab"
NAMD="namd2"

pp=./confs
FFCONV_PATH=`which ffconv.py | sed "s:/ffconv.py::g"`

lNAB='no'
lNAMD='yes'
#lgeometric='no'

# Currently we test the AMBER topology via NAMD, because: 
# 1) AMBER is not free
# 2) AmberTools (NAB) do not support many important fields of AMBER topology file

prog=AMBER

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
    --nab_cmd        - [default: $NAMD] executable for namd
    --namd_cmd       - [default: $NAB] executable for NAB
    --lnab           - [default: $lnab] perform energy evaluation with NAB of AmberTools ? yes/no
    --lnamd          - [default: $lnamd] perform energy evaluation with NAMD ? yes/no
    --scee           - custom scale constant for 1-4 electrostatic interactions 
    --scnb           - custom scale constant for 1-4 Lennard-Jones interactions
    --lgeometric     - [default: get from mol.top] use geometric combining rules ? yes/no
    --ambercrd       - file with coordinates in Ambedr format
    --lskipcheck     - skip energy evaluation
    -h               - print help
"

HelpMessage="

HELP MESSAGE FOR $0

The script performs a simple energy calculation for Amber topology file with NAMD or NAB. 

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
    --nab_cmd)
      NAB=$2; shift 2 
    ;;
    --namd_cmd)
      NAMD=$2; shift 2 
    ;;
    --lnab)
      lNAB=$2; shift 2 
    ;;
    --lnamd)
      lNAMD=$2; shift 2 
    ;;
    --scee)
      scee=$2; shift 2 
    ;;
    --scnb)
      scnb=$2; shift 2 
    ;;
    --lgeometric)
      lgeometric=$2; shift 2 
    ;;
    --ambercrd)
      ambercrd=$2; shift 2 
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
      echo >&2 -e "--- !Error in [ "$0" ]: the command line argument [$1] not known. Aborting.\n" ; exit 1;
    ;;
  esac
done


###------------------------------------------------------------------------------
### Checks

echo >&2 -e "--- !Note in ["$0"]: at the moment the energy for a Amber topology file can be calculated with NAMD package or NAB program of AmberTools package. Amber itself is not used because it is not free for academia. The Amber topology file is NOT fully recognized by NAMD or NAB.
\n---NAB:\n
1) ignores completely the SCEE and SCNB parameters of the topology file (and does not allow to set them explicitly).\n
2) Chamber type of the Amber topology is not supported as well.\n
Thus, so far only AMBER Force Field can be checked with NAB.
\n---NAMD:\n
1) ignores the SCEE and SCNB parameters of the topology file but allows to set them explicitly.\n
So far AMBER Force Field can be checked with NAMD fully, OPLSAA and CHARMM - only nonbonded interactions (including 1-4 pairs).\n
\n
Therefore we prefer to use NAMD for energy comparison. 
Continue.\n" ; 

if [ "$lNAMD" == 'yes' -a "$lNAB" == 'yes' ]; then
    echo >&2 -e "--- !Error in ["$0"]: lNAB and lNAMD must not be both yes at the same time. Aborting.\n" ; exit 1;
fi

if [ "$lNAMD" != 'yes' -a "$lNAB" != 'yes' ]; then
    echo >&2 -e "--- !Error in ["$0"]: lNAB and lNAMD must not be both no at the same time. Aborting.\n" ; exit 1;
fi


if [ "$lNAMD" == "yes" ]; then
    prog="NAMD"
fi

if [ "$lNAB" == "yes" ]; then
    prog="NAB"
fi


templdir=$FFCONV_PATH/templates/$prog
scriptsdir=$FFCONV_PATH/scripts


if [ "$lskipcheck" != 'yes' ]; then

mkdir $d

for i in `ls $pp/*.pdb`; do
    conf=`python -c "t=\"$i\".split('/'); t=t[-1]; t=t.split('.'); t=t[0]; print t"`

    mkdir $d/$conf 
    p=`pwd`; cd $d/$conf

     if [ "$lNAB" == 'yes' ]; then
        ftop=../../${mol}.top
        fpdb=$i

        FillVariablesToTemplate.sh FONAME=singlepoint.nab TEMPLATE=${templdir}/template_singlepoint.nab TMPPDB=${fpdb} TMPTOP=${ftop}
        
        cmd="$NAB -o nab.x  singlepoint.nab &> NABOUT"
        eval $cmd ; if [ $? -ne 0 ]; then { echo >&2 -e "--- !Error: [$cmd] failed. Aborting.\n" ; exit 1; } fi
        
        cmd="./nab.x > nab.out"
        eval $cmd ; if [ $? -ne 0 ]; then { echo >&2 -e "--- !Error: [$cmd] failed. Aborting.\n" ; exit 1; } fi
        
        $scriptsdir/get_${prog}_energy.sh nab.out > E.dat
     fi

     if [ "$lNAMD" == 'yes' ]; then
        ftop=../../${mol}.top
        fpdb=$i

        if [ -z "$scee" ]; then
            # scee not given as cmd argument => try to get it from top file
            scee=`cat ${ftop} | awk '/%FLAG SCEE_SCALE_FACTOR/ {flag=1;next} /%FLAG/{flag=0} flag {print}' | tail -n +2 | tail -n 1 | awk '{print $1}'`
        fi
        if [ -z "$scee" ]; then
            # scee stil empty => error
            echo >&2 -e "--- !Error in ["$0"]: scee not given as cmd argunment and not present in topology file [at least one value]. Aborting.\n" ; exit 1;
        fi
        fudgeQQ=`python -c "print 1.0/$scee"`

        if [ -z "$scnb" ]; then
            # scee not given as cmd argument => try to get it from top file
            scnb=`cat ${ftop} | awk '/%FLAG SCNB_SCALE_FACTOR/ {flag=1;next} /%FLAG/{flag=0} flag {print}' | tail -n +2 | tail -n 1 | awk '{print $1}'`
        fi
        if [ -z "$scnb" ]; then
            # scee stil empty => error
            echo >&2 -e "--- !Error in ["$0"]: scnb not given as cmd argunment and not present in topology file [at least one value]. Aborting.\n" ; exit 1;
        fi
        
        if [ -z "$lgeometric" ]; then
            # lgeometric not given as cmd argument => try to get it from top file
            CombRule=`cat ${ftop} | grep CombRule | awk '{print $3}'`

            if [ "$CombRule" == 'Geometric' ]; then
                lgeometric='yes'
            else
                lgeometric='no'
            fi
        fi

        if [ ! -f "$conf".crd ]; then
            # Amber CRD file do not exist for given PDB
            # Convert PDB to CRD
            echo >&2 -e "--- !Note in ["$0"]: no ${conf}.crd present. Generating it from PDB. Continue.\n" ;
            $scriptsdir/pdb2ambercrd.py ${i} > ${conf}.crd
        fi
        
        # TODO
        FillVariablesToTemplate.sh FONAME=namd.conf  TEMPLATE=${templdir}/template_namd.conf TMPCRD=${conf}.crd TMPTOP=${ftop} TMPFUDGEQQ=${fudgeQQ} TMPSCNB=${scnb} TMPVDWGEOMETRICSIGMA=${lgeometric}
        
        cmd="$NAMD namd.conf > namd.out"
        eval $cmd ; if [ $? -ne 0 ]; then { echo >&2 -e "--- !Error: [$cmd] failed. Aborting.\n" ; exit 1; } fi
        
        $scriptsdir/get_${prog}_energy.sh namd.out > E.dat
     fi


    cd $p
done


fi # lskipcheck


### Files to store to DB

#if [ ! -f "${mol}.crd" ]; then
#    echo >&2 -e "--- !Note in ["$0"]: no ${mol}.crd present. Copy the last entry from the confpath. Continue.\n" ;
#    cp ${pp}/${conf}.crd ${mol}.crd
#fi
if [ ! -f "${mol}.pdb" ]; then
    echo >&2 -e "--- !Note in ["$0"]: no ${mol}.pdb present. Copy the last entry from the confpath. Continue.\n" ;
    cp ${pp}/${conf}.pdb ${mol}.pdb
fi
    

echo "${mol}.pdb
${mol}.top
" >> MOLFILES







