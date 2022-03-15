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

# This script parses charmm outp file and gets the energy terms.
# The input file is 'charmm.out'

import sys
from FactorToConvertUnits import *

fn='charmm.out'
f=open(fn,'r')

D={}
cnt=-1
fields=[]
for line in f:
    ln=line.split()
    if len(ln) < 2: continue
    if line[0:5] == 'ENER ' or line[0:5] == 'ENER>':
        if ':' in ln[1]:
            try: del ln[ln.index('Eval#')]
            except: pass

            for v in ln[2:]:
                D[v]=float('NaN')
                fields.append(v)
        
        if '>' in ln[1] or '>' in ln[0]:
            #if ln[0] == 'ENER>': del ln[1]
            #except: pass
            for v in ln[2:]:
                cnt+=1
                D[fields[cnt]]=float(v)

#ENER ENR:  Eval#     ENERgy      Delta-E         GRMS
#ENER INTERN:          BONDs       ANGLes       UREY-b    DIHEdrals    IMPRopers
#ENER EXTERN:        VDWaals         ELEC       HBONds          ASP         USER
# ----------       ---------    ---------    ---------    ---------    ---------
#ENER>        0    -13.94836      0.00000     23.90908
#ENER INTERN>        5.23220     44.76648      1.67958      9.66305      0.00113
#ENER EXTERN>        9.62902    -84.91981      0.00000      0.00000      0.00000
# ----------       ---------    ---------    ---------    ---------    ---------

D2=D

lst=['E_POT','E_BONDS','E_ANGLES','E_DIHEDRALS','E_IMPROPERS','E_VDW','E_ELEC']

D={}
fact=FactorToConvertUnits('kcal/mol','kJ/mol')

if 'BONDs' in D2.keys(): D['E_BONDS'] = D2['BONDs']*fact
if 'UREY-b' in D2.keys(): 
    E_UB = D2['UREY-b']
else:
    E_UB = 0.0
    
if 'ANGLes' in D2.keys(): D['E_ANGLES'] = (D2['ANGLes'] + D2['UREY-b'] )*fact

if 'DIHEdrals' in D2.keys(): D['E_DIHEDRALS'] = D2['DIHEdrals']*fact
if 'IMPRopers' in D2.keys(): D['E_IMPROPERS'] = D2['IMPRopers']*fact

D['E_VDW'] = D2['VDWaals']*fact
D['E_ELEC'] = D2['ELEC']*fact
D['E_POT'] = D2['ENERgy']*fact

#sys.stderr.write(str(D2))
for k,v in D.iteritems():
    print k, v




