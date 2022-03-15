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

# 28 oct 2014. Bug fix: lines 66,67. Correct construction of atom types. In ffld_server output atom types do not distiguish uniquely vdw parameters 

### TODO: 
# Documentation        


def ReadFFLDServerOutput(fn):
    # ffld_server changes the names of the atoms from the PDB file: such that a name has element characters and the indext of the atom in molecule
    #
    def TrimAtomNumbers(L):
            tmp=[]; 
            for val in L: tmp.append( int(re.sub('[^0-9]','',val)))
            return tmp

    try: 
        tmp=fn.split('/'); tmp=tmp[-1].split('.') ;MoleculeName=tmp[0]
    except:
        MoleculeName='none'
     
    Mol = Molecule(MoleculeName)
    Mol.FFName='oplsaa'
    Mol.GeneratedBy='ffld_server'
    UnitsIn={ 'Energy':'kcal', 'Distance':'angstr', 'Angle':'rad', 'Number':'mol' }
    Mol.Units = UnitsIn
  
    # Constants
    fudgeLJ = FFvars['OPLS']['fudgeLJ'] # OPLS
    
    ifile=open(fn,'r')
    while 1:
        line=ifile.readline()
        if not line: break
        ln=line.split()
        if len(ln)<5: continue
        if ln[0]=='atom' and ln[1]=='type' and ln[2]=='vdw':
            line=ifile.readline()
            iAtom = 0
            while 1:
                iAtom += 1
                line=ifile.readline()
                ln=line.split()
                if (len(ln)<3): break
                A = Atom(iAtom)
                A.Name = ln[0]
                attype = ln[1]
                #atsymbol = TrimElementName(ln[3])
                atsymbol = ln[2]
                A.Type = atsymbol +'_'+ attype
                A.Charge = float(ln[4])
                A.Func = 'LJ_sig_eps'
                A.p =  [ float(ln[5]) , float(ln[6]) ] 
                A.p_14 = [ float(ln[5]), float(ln[6])*fudgeLJ ] #Note!!! For OPLS we use fudgeLJ = 0.5 for 1-4 pairs. For FuncType LJ_sig_eps, we need just multiply eps by fudgeLJ. Note!!! if you use these 1-4 parameters set gen-pairs to no.
                ### [ angstr, kcal/mol ]               
                
                A.Quality = ln[7]
                A.Description = " ".join(ln[8:])
                #A.Comment = ""
                
                Mol.Atoms.append(copy(A))

                
        elif ln[0]=='Stretch' and ln[1]=='k' and ln[2]=='r0':
            iBond = 0
            while 1:
                iBond += 1
                line = ifile.readline()
                ln = line.split()
                if (len(ln)<3): break
                
                B = Bond(iBond)
                B.AtNames = [ ln[0], ln[1] ] 
                B.AtNumbers = TrimAtomNumbers(B.AtNames)
                B.Func = 'Harmonic'
                B.p = [ float(ln[3]), float(ln[2])*2.0 ] # Convert const: Kgmx = 2*Kffld
                ### [ angstr, kcal/mol/angstr^2 ]               
                B.Quality = ln[4]
                #B.Description = ""
                B.Description = " ".join(ln[8:])
                Mol.Bonds.append(copy(B))
 
        elif ln[0]=='Bending' and ln[1]=='k' and ln[2]=='theta0':
            iAngle = 0
            while 1:
                iAngle += 1
                line = ifile.readline()
                ln = line.split()
                if (len(ln)<3): break
                A = Angle(iAngle)
                A.AtNames = [ ln[0], ln[1], ln[2] ] 
                A.AtNumbers = TrimAtomNumbers(A.AtNames)
                A.Func = 'Harmonic'
            ### ! Mind here: the units of angles in written gromacs topologies are DEGREEs, while in the force constants RADs. 
            ### Internal units of gromacs (which used in the code) for angles are RADs.
                A.p = [ float(ln[4])*FactorToConvertUnits('degree','rad'), float(ln[3])*2.0 ] # Convert const: Kgmx = 2*Kffld
                ### [ rad , kcal/mol/rad^2 ] 
                A.Quality = ln[5]
                #A.Description = ""
                A.Description = " ".join(ln[7:])
                Mol.Angles.append(copy(A))

        elif ln[0]=='proper' and ln[1]=='Torsion' and ln[2]=='V1':
            iTor = 0
            while 1:
                iTor += 1
                line = ifile.readline()
                ln = line.split()
                if (len(ln)<3): break
                T = Torsion(iTor)
                T.AtNames = [ ln[0], ln[1], ln[2], ln[3] ] 
                T.AtNumbers = TrimAtomNumbers(T.AtNames)
                T.Func = 'Fourier'
                T.p = [ float(ln[4]), float(ln[5]), float(ln[6]), float(ln[7]) ] 
                ### [ 4 times kcal/mol ] But radians must be there also!              
                T.Quality = ln[8]
                #T.Description = ""
                T.Description = " ".join(ln[10:])
                Mol.Torsions.append(copy(T))
        
        elif ln[0]=='improper' and ln[1]=='Torsion' and ln[2]=='V2':
            iImp = 0
            while 1:
                iImp += 1
                line = ifile.readline()
                ln = line.split()
                if (len(ln)<3): break
                I = Improper(iImp)
                I.AtNames = [ ln[0], ln[1], ln[2], ln[3] ] 
                I.AtNumbers = TrimAtomNumbers(I.AtNames)
                I.Func = 'Fourier'
                I.p = [ 0.0, float(ln[4]), 0.0, 0.0 ]  
                ### [ 4 times kcal/mol ] But radians must be there also!              
                I.Quality = ln[5]
                #I.Comment = ""
                I.Description = " ".join(ln[6:])
                Mol.Impropers.append(copy(I))


    ifile.close()
    
    return Mol


