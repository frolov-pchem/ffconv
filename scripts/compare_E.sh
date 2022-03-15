#/bin/bash
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

# TODO:
#        required scripts checks

### Default values
pp=confs
lsuccess='yes'
FFCONV_PATH=`which ffconv.py | sed "s:/ffconv.py::g"`

KnownProgs=( "GROMACS" "CHARMM" "FFLDSERVER" "AMBER")


BriefUsage="--------------------------------------------------------------------------
BRIEF USAGE (for $0)
    \$FFCONVPATH/scripts/compare_E.sh --ffconvpath \$FFCONVPATH --fromprog FFLDSERVER --toprog GROMACS 

"

HelpMessage="
-----------------------------------------------------------------------------
Help message for $0

The script performs comparison of the energy terms with the help of compare_E.py script. It is used as a subroutine of check_conversion.sh script. For this to work corrrectly you should have already performed energy evaluations for each PDB configuration and have the files E.dat with energy contribution energies which will be compared for the two programs. If each energy term difference is less than $tr kJ/mol then the FF files conversion is considered to be successful: line with SUCCESS is writen to standard output.

${BriefUsage}

OPTIONS
    --confpath     - path to PDB files with different conformations of the molecule under consideration 
    --ffconvpath   - [default: ${FFCONV_PATH} - deduced from location of ffconv.py] path to ffconv (and the utilities)
    --fromprog     - program to convert from
    --toprog       - program to convert to
    --compareonly  - [default: all] specify energy terms to compare given separated by comma. E.g. E_POT,E_VDW
    -h             - print help and exit


KNOWN PROGS
${KnownProgs[*]}

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
    --ffconvpath)
      FFCONV_PATH=$2; shift 2 
    ;;
    --fromprog)
      fromprog=$2; shift 2 
    ;;
    --toprog)
      toprog=$2; shift 2 
    ;;
    --compareonly)
      compareonly=$2; shift 2 
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
### Parse command line arguments
###------------------------------------------------------------------------------


### Vars
scriptsdir=${FFCONV_PATH}/scripts
d1=$fromprog
d2=$toprog



### Checks

if [ ! -f "`find $pp/*.pdb`" ]; then
        echo >&2 -e "--- !Error in $0: no PDB files present in path [${pp}]. Expect to have multiple pdb files with different conformations of the stadied molecule. Nothing to do. Exit." ; exit 1
fi

if [ ! -f "$scriptsdir/compare_E.py" ]; then
        echo >&2 -e "--- !Error in $0: no ${scriptsdir}/compare_E.py found. Can not perform energy comparison. Exit." ; exit 1
fi

### Cycle over all availbale conformation PDB present. Should be the same as the conformers for which energies calculated.
for i in `ls $pp/*.pdb`; do
    conf=`python -c "t=\"$i\".split('/'); t=t[-1]; t=t.split('.'); t=t[0]; print t"`

    if [ ! -z "$compareonly" ]; then
        cokey="--compareonly ${compareonly}"
    else
        cokey=""
    fi

    # Perform energy comparison
    leq=`$scriptsdir/compare_E.py --fn1 $d1/$conf/E.dat --fn2 $d2/$conf/E.dat  ${cokey}   2> ${conf}_compare_E.log`
    echo ${leq}
    
    if [ $? -ne 0 ]; then { echo >&2 -e "--- !Error: [compare_E.py] failed. Aborting.\n" ; exit 1; } fi
    
    if [ "$leq" != "yes" ]; then
        echo >&2 -e "--- !Warning in $0: energies for $conf do not match. See ${conf}_compare_E.log for details. Continue."
        lsuccess='no'
    fi
done


### Print the result
if [ "$lsuccess" != "yes" ]; then
    echo -e "\n--- !Warning in $0: FF files conversion from $d1 format to $d2 format is not correct. Think what to do...\n"
else
    echo -e "\n--- SUCCESS in $0: FF files conversion from $d1 format to $d2 format correct.\n"
fi

