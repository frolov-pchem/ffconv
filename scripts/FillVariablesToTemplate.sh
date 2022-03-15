#!/bin/bash
#
# Written:       by Andrey Frolov     2011, MIS MPI, Leipzig, Germany
# Last modified: by Andrey Frolov apr 2013, ISC RAS, Ivanovo, Russia
#
#   Copyright 2013 Andrey I. Frolov
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


UsageMSG="
USAGE
$0 FONAME=fileout.dat TEMPLATE=template_for_fileout.dat TMPVAR1=some_value1 TMPVAR2=some_value2 ... TMPVARN=some_valueN
(e.g. $0 FONAME=plumed.dat TEMPLATE=template_plumed.dat TMPCV=10.0 TMPKS=1000.0)

The script $0 copies TEMPLATE file to FONAME file and substitutes all matches of TMPVAR1 by some_value1, all matches of TMPVAR2 by some_value2, etc.
FONAME and TEMPLATE are abligatory to be specified.
"

# Check if cmd argumets given
if [ -z "$*" ]; then
    echo >&2 -e "$UsageMSG"
    echo >&2 -e "!Error in $0: no cmd line arguments given. Aborting.\n"; exit 1;
fi

# Check if Help message requested
if [ "$1" == "-h" ]; then
    echo -e "$UsageMSG"
    exit 0;
fi



###------------
### Get FONAME and TEMPLATE files names
for i in "$@" ; do
        var="`echo $i | awk -F '=' '{print $1}'`"
        val="`echo $i | awk -F '=' '{print $2}'`"
        if [ "$var" == "FONAME" ]; then
                FONAME=$val
        elif [ "$var" == "TEMPLATE" ]; then
                TEMPLATE=$val
        fi
done

###-----------
### Checks

if [ -z "$FONAME" ]; then
    echo >&2 -e "!Error in $0: FONAME key is not set. Aborting. "; exit 1;
fi

if [ -z "$TEMPLATE" ]; then
    echo >&2 -e "!Error in $0: TEMPLATE key is not set. Aborting. "; exit 1;
fi

if [ ! -f "$TEMPLATE" ]; then 
    echo "!Error in $0: template file [$TEMPLATE] doesn't exist. Aborting."; exit 1;
fi

if [ -f "$FONAME" ]; then 
    echo "!Warning in $0: output file [$FONAME] already exist. Make a backup and overwrite.";
    # Check if bk.sh is in the PATH 
    cmd=bk.sh
    command -v "$cmd" >/dev/null 2>&1 || { echo -e >&2 "!Error in $0: script/program [$cmd] isn't in the PATH, but you requested backuping a file. Aborting.\n"; exit 1; }
    bk.sh "$FONAME"
fi


###-----------
### Copy template to output file
cp $TEMPLATE $FONAME


###-----------
### Substitute temporary fields in file
for i in "$@" ; do
	var="`echo $i | awk -F '=' '{print $1}'`"
	val="`echo $i | awk -F '=' '{print $2}'`"
	if [ "$var" != "FONAME" -a "$var" != "TEMPLATE" ]; then

                if [ -z "$var" ]; then
                    echo >&2 -e "!Error in $0: a template field is empty, which was requested to be substituted by [$val]. Aborting."; exit 1
                fi
                
                if [ -z "$val" ]; then
                    echo >&2 -e "!Warning in $0: a value for template field [$var] is empty. Continue.";
                fi

                Nl=`cat $FONAME | grep "$var" | wc -l`
                if [ "$Nl" -gt "0" ]; then
		    #echo sed -i "s@${var}@${val}@g" $FONAME
		    sed -i "s+${var}+${val}+g" $FONAME
                else
                    echo >&2 -e "!Warning in $0: template file [$TEMPLATE] has no fields [${var}], whcih were requested to be substituted by [$val]. Continue.";
                fi
	fi
done




