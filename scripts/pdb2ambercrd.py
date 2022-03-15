#!/usr/bin/env python
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

# The script converts pdb to amber coor file.
# Input file is taken as a cmdline argument.

# Andrey I. Frolov, Jun 2013, ISC RAS, Ivanovo

import sys, numpy

def PrintAmberTopFieldValues(lst,step,fmt):
    i=0
    s=''
    while i < len(lst):
        step = min(step,len(lst)-i)
        fmt2 = step * fmt+"\n"
        s += fmt2 % tuple(lst[i:i+step])  
        i+=step
    return s


fname=sys.argv[1]
ifile=open(fname,"r")
lst=[]
for line in ifile:
    ln=line.split()
    if len(ln)==0: continue
    if line[0:4]=='ATOM' or line[0:4]=='HETATM':
        lst +=  [ float(line[30:38]), float(line[38:47]), float(line[47:56]) ]

s = 'Converted from '+fname+' by '+sys.argv[0]+'\n'
s += PrintAmberTopFieldValues([len(lst)/3],1,"%5i")
s += PrintAmberTopFieldValues(lst,6,"%12.7f")

print s

#ATOM      1  O1  PARA    1       3.537   1.423  -0.000  1.00  0.00
#ATOM      2  C1  PARA    1       3.620   2.617  -0.047  1.00  0.00



