#!/usr/bin/env python
import sys
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

# Last modified: 22 may 2014

# The script performs comparison of energy terms stored in two different files

# Threshold of the energy difference (to account for minor differences between different codes)
tr=0.1 #kJ/mol

###------------------------------------------------------
### Options parser
from optparse import OptionParser

usage = "usage: %prog [options] arg1 arg2"
parser = OptionParser(usage=usage)
parser.add_option("--fn1", 
                  default='',
                  dest="fn1", type="string", 
                  help="input file name 1")
parser.add_option("--fn2", 
                  default='',
                  dest="fn2", type="string", 
                  help="input file name 2")
parser.add_option("--compareonly", 
                  default='',
                  dest="compareonly", type="string", 
                  help="test only listed by comma energy terms")

(opts, args) = parser.parse_args()

### Set options
fn1 = opts.fn1
fn2 = opts.fn2
Elist = opts.compareonly.split(',')


def ReadFile(fn):
    f=open(fn,'r')
    D={}
    for line in f:
        ln=line.split()
        D[ln[0]]=float(ln[1])
    return D



D1 = ReadFile(fn1)
D2 = ReadFile(fn2)

if len(D1.keys()) > len(D2.keys()):
    DL = D1
    DS = D2
    fnL = fn1
    fnS = fn2
else:
    DL = D2
    DS = D1
    fnL = fn2
    fnS = fn1
    

leq='yes'

if Elist == ['']:

    if not 'E_POT' in DS.keys():
        sys.stderr.write('--- !Error: the key [E_POT] must be present in the energy file. Could not find it in ['+fnS+']. Aborting.\n')
        leq='no'
        print leq
        sys.exit(1)

    if not 'E_POT' in DL.keys():
        sys.stderr.write('--- !Error: the key [E_POT] must be present in the energy file. Could not find it in ['+fnL+']. Aborting.\n')
        leq='no'
        print leq
        sys.exit(1)

else:

    for v in Elist:
        if not v in DS.keys():
            sys.stderr.write('--- !Error: the key ['+v+'] must be present in the energy file. Could not find it in ['+fnS+']. Aborting.\n')
            leq='no'
            print leq
            sys.exit(1)
        if not v in DL.keys():
            sys.stderr.write('--- !Error: the key ['+v+'] must be present in the energy file. Could not find it in ['+fnL+']. Aborting.\n')
            leq='no'
            print leq
            sys.exit(1)


used_keys=[]
for k,v in DL.iteritems():
    
    if Elist != ['']:
        if not k in Elist:
            sys.stderr.write('--- !Warning: the key ['+k+'] from file ['+fnL+'] not used in energy comparison. Continue.\n')
            continue
            
    if not k in DS.keys():
        sys.stderr.write('--- !Warning: the key ['+k+'] from file ['+fnL+'] not used in energy comparison. Continue.\n')
        continue
    else:
        sys.stderr.write('the key ['+k+'] from file ['+fnL+'] used in energy comparison.\n')
        
        
    used_keys += [ k ] 
    if abs(DS[k] - DL[k]) > tr:
        leq='no'
        sys.stderr.write('diff   '+k+'   '+str(abs(DS[k]-DL[k])))
        sys.stderr.write('--- !Warning: the differencers in energies of ['+k+'] exceed treshold (kJ/mol) ['+str(tr)+']. Continue.\n')



for k,v in DS.iteritems():
    if not k in used_keys:
        sys.stderr.write('--- !Warning: the key ['+k+'] from file ['+fnS+'] not used in energy comparison. Continue.\n')
    else:
        sys.stderr.write('the key ['+k+'] from file ['+fnS+'] used in energy comparison.\n')
    
# Print the result of comparison: are identical ? 
print leq       


