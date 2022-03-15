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

f=$1
suf=$(date "+%F_%T") 
suf=`echo $suf | sed "s/\:/\./g"`

fout=${f}.bk.${suf}

### Check
if [ -z "$f" ]; then
    echo >&2 -e "
The script makes a backup copy of a file or a directory (of cause recursively!).
USAGE: bk.sh file
USAGE: bk.sh dir 

OUTPUT
file.bk.date
dir.bk.date
"
    echo >&2 -e "!Error in $0: no cmd argument given. Aborting."; exit 1
fi


if [ -f "$f" ]; then
    echo "---Back up file"
    cp ${f} ${fout}
    if [ -f  "${fout}" ]; then 
        echo "---File: ${fout} created"
    else
        echo "---!!!Error: backup file ${fout} NOT created by some reason. Aborting."; exit 1
    fi
elif [ -d "$f" ]; then
    echo "---Back up directory"
    cp -r ${f} ${fout}
    if [ -d  "${fout}" ]; then 
        echo "---Directory: ${fout} created"
    else
        echo "---!!!Error: backup directory ${fout} NOT created by some reason. Aborting."; exit 1
    fi
else
    echo "!Error in $0. The file "$f" is not a file or a directory. Aborting."; exit 1
fi

exit 0

