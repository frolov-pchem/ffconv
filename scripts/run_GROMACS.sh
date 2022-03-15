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
GROMPP="gmx grompp"
MDRUN="gmx mdrun"
GENERGY="gmx energy"

pp=./confs
FFCONV_PATH=`which ffconv.py | sed "s:/ffconv.py::g"`


d=GROMACS

BriefUsage="--------------------------------------------------------------------------
BRIEF USAGE (for $0)

OPTIONS

    --mol             - name of the molecule (defines the names of input files)
    --confpath        - [default: ${confpath}] path to where PDB file(s) with coordinates of molecule located.
                        if many present and lskipcheck != yes: then energies of configuration from all pdb files 
                        are used to check the correctness of conversion.
    --ffconvpath      - [default: \$FFCONVPATH] path to ffconv distribution. If FFCONV not defined then the path set 
                        to location of ffconv.py executable.
    --gmx-ff-itp      - gromacs ITP force field file - common entries for all molecules 
    --gmx-ff-defaults - [default: defaults.itp] gromacs ITP force field file - defaults section  
    --grompp_cmd      - [default: ${GROMPP}] grompp executable
    --mdrun_cmd       - [default: ${MDRUN}] mdrun executable
    --genergy_cmd     - [default: ${GENERGY}] g_energy executable
    --lskipcheck      - skip energy evaluation
    -h                - print help

"

HelpMessage="-------------------------------------------------------------------------
Help message for $0

DESCRIPTION

This script performs energy calculation for single molecules configurations with Gromacs

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
    --grompp_cmd)
      GROMPP=$2; shift 2 
    ;;
    --mdrun_cmd)
      MDRUN=$2; shift 2 
    ;;
    --genergy_cmd)
      GENERGY=$2; shift 2 
    ;;
    --gmx-ff-defaults)
      gmxffdefaults=$2; shift 2 
    ;;
    --gmx-ff-itp)
      gmxffitp=$2; shift 2 
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
      echo >&2 -e "!Error in [ "$0" ]: the command line argument [$1] not known. Aborting.\n" ; exit 1;
    ;;
  esac
done


### Set default 
    if [ -z "$gmxffdefaults" ]; then
        #cp ../defaults.itp .
        gmxffdefaults=`pwd`/defaults.itp
    fi    



### Skip energy evaluation
if [ "$lskipcheck" != "yes" ]; then

templdir=$FFCONV_PATH/templates/GROMACS
scriptsdir=$FFCONV_PATH/scripts


mkdir $d
p=`pwd`; cd $d

    
    ### Gen Topology
s1="
;;; LOAD DEFAULTS
#include \"${gmxffdefaults}\" 

"

if [ ! -z "$gmxffitp"  ]; then
s2="

;;; LOAD FF *.itp
#include \"${gmxffitp}\" 

"
fi

s3="
;;; LOAD MOLECULES *.itp
#include \"${p}/${mol}_AT.itp\" 

;;; LOAD MOLECULES *.itp
#include \"${p}/${mol}.itp\" 

[ system ]                             
; Name                                 
System: ${mol}

[ molecules ]
 ${mol}  1
"

    echo "$s1 $s2 $s3" > topol.top

#    FillVariablesToTemplate.sh FONAME=topol.top TEMPLATE=${templdir}/template_singlemol_topol.top TMPMOL=${mol} TMPFFITP=${gmxffitp} TMPDEFAULTS=${gmxffdefaults} 
    

    run=em
    cp ${templdir}/${run}.mdp .
    

cd $p


### Loop over all single molecule configurations
for i in `ls $pp/*.pdb`; do
    conf=`python -c "t=\"$i\".split('/'); t=t[-1]; t=t.split('.'); t=t[0]; print t"`

    mkdir $d/$conf 
    p=`pwd`; cd $d/$conf

        # Energy minimization
        
        name=em
        mdp=../${run}.mdp
        itop=../topol.top
        #igro=$pp/$i
        igro=$i
        cmd="$GROMPP -f ${mdp} -p ${itop} -c ${igro} -o ${name} -maxwarn 6 &> GROMPPOUT_${name}"
        eval $cmd ; if [ $? -ne 0 ]; then { echo >&2 -e "--- !Error: [$cmd] failed. Aborting.\n" ; exit 1; } fi
        cmd="$MDRUN -deffnm ${name}  &> MDRUNOUT_${name}"
        eval $cmd ; if [ $? -ne 0 ]; then { echo >&2 -e "--- !Error: [$cmd] failed. Aborting.\n" ; exit 1; } fi
        

# TODO list of fileds        
        $GENERGY -f em.edr -b 0 -e 0 -o em.xvg > energy.out 2>/dev/null <<EOL
Bond
U-B
Angle
Proper-Dih.
Fourier-Dih.
Ryckaert-Bell.
Improper-Dih.
LJ-14
Coulomb-14
LJ-(SR)
Coulomb-(SR)
Potential
EOL

        if [ $? -ne 0 ]; then { echo >&2 -e "--- !Error: [$GENERGY] failed. Aborting.\n" ; exit 1; } fi

        $scriptsdir/get_GROMACS_energy.py energy.out > E.dat

    cd $p
done


fi # lskipcheck


### Files to store to DB

if [ ! -f "${mol}.pdb" ]; then
    echo >&2 -e "--- !Note in ["$0"]: no ${mol}.pdb present. Copy the last entry from the confpath. Continue.\n" ;
    cp ${pp}/${conf}.pdb ${mol}.pdb
fi

echo "${mol}.pdb
${mol}_AT.itp
${mol}.itp
" >> MOLFILES

echo "$gmxffdefaults
" >> FFFILES
if [ ! -z "$gmxffitp"  ]; then
echo "$gmxffitp
" >> FFFILES
fi

