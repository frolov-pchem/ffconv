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

# This script parses the gromacs g_energy output file and greps the energy terms.
# The input file is taken as cmdline argument.

# TODO
# 1) documentation
# 2) make code optimization ?
# 3) Missing keys in dictionary 
# 4) factor to convert units ? get from unitconv.py

import sys

list=['Bond','U-B','Angle','Proper-Dih.','Fourier-Dih.','Improper-Dih.','LJ-14','Coulomb-14','LJ-(SR)','Coulomb-(SR)','Potential','Ryckaert-Bell.']
lst=['E_POT','E_BONDS','E_ANGLES','E_DIHEDRALS','E_IMPROPERS','E_VDW','E_ELEC']

fn=sys.argv[1]

f=open(fn,'r')

D={}
for line in f:
    ln=line.split()
    if len(ln) < 2: continue
    if ln[0]+'-'+ln[1] in list:
        D[ln[0]+'-'+ln[1]]=float(ln[2])
        continue
    if ln[0] in list:
        D[ln[0]]=float(ln[1])


D2=D

D={}

# kJ/mol are default values.
fact=1.0

D['E_POT'] = D2['Potential']*fact
D['E_BONDS'] = D2['Bond']*fact
try:
    D['E_ANGLES'] = D2['U-B']*fact
except:
    D['E_ANGLES'] = D2['Angle']*fact

if 'Proper-Dih.' in D2.keys(): 
#    D['E_TORSIONS'] = D2['Proper-Dih.']*fact
    D['E_DIHEDRALS'] = D2['Proper-Dih.']*fact

if 'Fourier-Dih.' in D2.keys(): 
#    D['E_TORSIONS'] = D2['Fourier-Dih.']*fact
    D['E_DIHEDRALS'] = D2['Fourier-Dih.']*fact

if 'Ryckaert-Bell.' in D2.keys(): 
#    D['E_TORSIONS'] = D2['Ryckaert-Bell.']*fact
    D['E_DIHEDRALS'] = D2['Ryckaert-Bell.']*fact

if 'Improper-Dih.' in D2.keys(): 
#    D['E_IMPROPERS'] = D2['Improper-Dih.']*fact
    D['E_DIHEDRALS'] += D2['Improper-Dih.']*fact

D['E_VDW'] = (D2['LJ-(SR)'] +  D2['LJ-14']  )*fact
D['E_ELEC'] = (D2['Coulomb-(SR)'] +  D2['Coulomb-14'] )*fact

D['E_BONDED'] = (D['E_DIHEDRALS'] +  D['E_ANGLES'] + D['E_BONDS'] )*fact



for k,v in D.iteritems():
    print k, v








