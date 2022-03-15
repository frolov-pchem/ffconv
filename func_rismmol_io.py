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
from molecule_class import *
from copy import copy
import re
#NaN=float('NaN')

### TODO: 
# Documentation        



#----------------------------------------------------------------------------------------------------------------------
def WriteRISMMOLMoleculeTopFile(FName,Mol,**kwargs):
#
#   EXAMPLE USAGE
#
#   DESCRIPTION
#   
#   INPUT
# Fname - name of the output file (string)
# Mol - a Molecule class object
#
#   REQUIREMENTS
#
#
    ### TODO: 
    
    NBFuncType = 'LJ_sig_eps'


    ###
    ### Necessary input checks
    ###
    if Mol.__class__.__name__ != 'Molecule':
            sys.stderr.write('--- !Error in function '+inspect.stack()[0][3]+': the Second function argument must be of class Molecule. Class of the given object: ['+Mol.__class__.__name__+']. Aborting.\n')
            sys.exit(1)

    if FName.__class__.__name__ != 'str':
            sys.stderr.write('--- !Error in function '+inspect.stack()[0][3]+': the First function argument must be of class str. Class of the given object: ['+FName.__class__.__name__+']. Aborting.\n')
            sys.exit(1)
  

    ###
    ### Checks
    ###
    if len(Mol.Atoms) < 1:
            sys.stderr.write('--- !Error in function '+inspect.stack()[0][3]+': the number of atoms for the molecule named ['+Mol.Name+'] = 0. Not allowed. Aborting.\n')
            sys.exit(1)
        

    ###
    ### Convert NB params: e.g from LJ_rmin2_eps tp LJ_sig_eps
    ###
    ConvertNonBondedParamsTo(Mol,NBFuncType)


    # Must be global!
    #MolFields = ['Atoms','Pairs','Bonds','Angles','Torsions','Impropers']
    
    ### AtomNBFuncType
    D_FT={}
    lst=getattr(Mol,'Atoms')
    D_FT['Atoms']=[]
    for F in lst:
        if not F.Func == 'none':
            D_FT['Atoms'].append(F.Func)
        else:
            sys.stderr.write('--- !Error in function '+inspect.stack()[0][3]+': in the molecule instance named ['+Mol.Name+'] the functional form for an atoms is none. The atom field is: ['+str(F.__dict__)+']. Aborting\n')
            sys.exit(1)

    # Mind! For RISMOL topology all atom types must have the same functional form LJ_sig_eps
    AtomNBFuncType = CheckIfAllElementsInArrayAreTheSameAndReturn(D_FT['Atoms'])
    if bool(AtomNBFuncType) == 0:
            sys.stderr.write('--- !Error in function '+inspect.stack()[0][3]+': in the Mol instance named ['+Mol.Name+'] the functional forms are not identical for some atoms! This is required by RISMOL topology file format. Aborting\n')
            sys.exit(1)
        

    ###
    ### Perform conversions
    ###
    
    # Convert units to GMX
    #D={'Energy':'kj','Angle':'rad','Distance':'angstr','Number':'mol'}
    ConvertUnitsTo(Mol,'toWHAM')

    ###
    ### Writing the ITP file
    ###
    of=open(FName,'w')
    string = ''

    for A in Mol.Atoms:

        if isnan(A.Charge): 
            sys.stderr.write('--- !Error in function '+inspect.stack()[0][3]+': the molecule instance named ['+Mol.Name+'] has atoms with charge set to NaN. The Atom instance ['+str(F.__dict__)+']. Aborting.\n ')
            sys.exit(1)

        if isnan(A.p[0]): 
            sys.stderr.write('--- !Error in function '+inspect.stack()[0][3]+': the molecule instance named ['+Mol.Name+'] has atoms with sigma set to NaN. The Atom instance ['+str(F.__dict__)+']. Aborting.\n ')
            sys.exit(1)

        if isnan(A.p[0]): 
            sys.stderr.write('--- !Error in function '+inspect.stack()[0][3]+': the molecule instance named ['+Mol.Name+'] has atoms with epsilon set to NaN. The Atom instance ['+str(F.__dict__)+']. Aborting.\n ')
            sys.exit(1)


        s=" %14.4f %14.4f %14.4f %14.8f %14.8f %14.8f " % (A.X, A.Y, A.Z, A.p[0], A.p[1], A.Charge) 
        s+='\n'
        string+=s
    string+='\n \n'

    of.write(string)
    return

