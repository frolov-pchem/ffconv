from molecule_class import *
#from copy import copy
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

NaN=float('NaN')

def ReadCharmmForceFieldFile(fn):
# TODO: get fudgeLJ from FF_global_dict ?
# TODO: get FFName from file or cmdline arg

    FF = ForceField('cgenff2b7')
    FF.Units = { 'Energy':'kcal', 'Distance':'angstr', 'Angle':'rad', 'Number':'mol' }
    FF.Family = 'CHARMM'
    FF.NREXCL = 3
    FF.fudgeLJ = 1.0
    FF.fudgeQQ = 1.0
    FF.GeneratedBy = 'paramchem'
    FF.CombRule = 'LorentzBerthelot'

    FieldsList=['ATOMS','BONDS','ANGLES','DIHEDRALS','IMPROPERS','NONBONDED','HBOND','END']

    ifile=open(fn,'r')
    field='none'
    D={}
    while 1:
        line=ifile.readline()
        if not line: break
        line=line[:-1]
        
        
        if '!' in line: 
            ln=line.split('!')
            line=ln[0]
            comment=ln[1]
        ln=line.split()
        if len(ln)==0: continue
        
        
	if ln[0][0]=='!' or ln[0][0]=='*': 
                continue
        
        if line[-1]=='-': 
            line+=ifile.readline()
            line=line[:-1]
        
        if ln[0] in FieldsList:
            field=ln[0]
        


#ATOMS
#!hydrogens
#MASS   256 HGA1     1.00800  ! alphatic proton, CH
#MASS   257 HGA2     1.00800  ! alphatic proton, CH2
        ### Mind! Assume that all MASS instances occur earlier than the NONBONDED section
        if field=='ATOMS':
            if ln[0]!='ATOMS':
                D[ln[2]] = {'Mass':float(ln[3]), 'Comment':comment } 


        if field=='NONBONDED':
            if ln[0]!='NONBONDED':
                A = AtomType()
                A.Type = ln[0]
                A.Mass = D[A.Type]['Mass']
                try: A.PeriodicTableNum = Elements[TrimElementName(A.Type)]['PeriodicTableNum']
                except: pass
                A.Charge = 0.0
                A.Func = 'LJ_Rmin2_eps'
                A.p = [ float(ln[3]),abs(float(ln[2])) ]
                try:
                    A.p_14 = [ float(ln[6]),abs(float(ln[5])) ]  
                except:
                    A.p_14 = copy(A.p) 
                #ln[2] # Rmin/2 [Angstr]
                #ln[3] # Eps [kcal/mol]
                #ln[4] # Rmin/2 [Angstr] for 1-4 interacrtions
                #ln[5] # Eps [kcal/mol] for 1-4 interactions
                
                A.Description = comment + ". "+D[A.Type]['Comment']
                FF.AtomTypes.append(copy(A))
       
        elif field=='BONDS':
            if ln[0]!='BONDS':
                if (len(ln)<3): break
                B = BondType()
                B.AtTypes = [ ln[0],ln[1] ] 
                B.Func = 'Harmonic'
                B.p = [ float(ln[3]), 2.0*float(ln[2]) ]  
                #ln[2] Kb, [kcal/mol/A^2]. Mind we need to convert to Gromacs convention: Kgmx = 2*Kcharmm
                #ln[3] r0, [A]
                B.Description = comment
                FF.BondTypes.append(copy(B))

        elif field=='ANGLES':
            if ln[0]!='ANGLES':
                A = AngleType()
                A.AtTypes = [ln[0],ln[1],ln[2]]
                ### ! Mind here: the units of angles in written gromacs topologies are DEGREEs, while in the force constants RADs. 
                ### Internal units of gromacs (which used in the code) for angles are RADs.
                A.p = [ float(ln[4])*FactorToConvertUnits('degree','rad'), float(ln[3])*2.0 ] # Convert const: Kgmx = 2*Kcharmm
                ### [ rad , kcal/mol/rad^2 ] 
                #ln[3] k_phi, [kcal/mol/rad^2]. Mind we need to convert to Gromacs convention: Kgmx = 2*Kcharmm
                #ln[4] phi0, [degrees]
                if len(ln) > 5:
                    try:
                        ln[5]=float(ln[5])
                        ln[6]=float(ln[6])
                    except:
                        ln[5]=0.0
                        ln[6]=0.0
                    
                    A.p.append( float(ln[6]) ) 
                    A.p.append( 2.0*float(ln[5]) )
                    #ln[5] k_ub, [kcal/mol/A^2]. Mind we need to convert to Gromacs convention: Kgmx = 2*Kcharmm
                    #ln[6] r_ub, [A]
                else: 
                    A.p.append( 0.0 ) # If no UB but harmonic set the parameters to zero to use UB potential form
                    A.p.append( 0.0 ) # If no UB but harmonic set the parameters to zero to use UB potential form
                A.Func = 'UreyBradley'
                A.Description = comment
                FF.AngleTypes.append(copy(A))

        elif field=='DIHEDRALS':
            if ln[0]!='DIHEDRALS':
                T = TorsionType()
                T.AtTypes = [ ln[0],ln[1],ln[2],ln[3] ] 
                T.Func = 'PeriodicMultiple'
         # #####pTn.append( [ float(ln[6])+180.0,float(ln[4]),float(ln[5]) ] ) 
                T.p = [ float(ln[6])*FactorToConvertUnits('degree','rad'),float(ln[4]),float(ln[5]) ]  
                #ln[4] K_phi, No need convert, the same functional form
                #ln[5] multiplicity, 
                #ln[6] phi0 - phase. Mind we need NOT (?) to convert to Gromacs convention: cis configuration = 0 angle, in Charmm it is cis too.
                #       ! Mind: need to converto to rads for consistency!
                T.Description = comment
                FF.TorsionTypes.append(copy(T))

        elif field=='IMPROPERS':
            if ln[0]!='IMPROPERS':
                I = ImproperType()
                I.AtTypes = [ ln[0],ln[1],ln[2],ln[3] ] 
                I.Func = 'Harmonic'
                #pIm.append( [ float(ln[6]) , 2.0*float(ln[4]) , float(ln[5]) ] ) 
                I.p = [ float(ln[6])*FactorToConvertUnits('degree','rad') , 2.0*float(ln[4]) ] 
                #ln[4] K_phi, 
                #ln[5] ignored (0), 
                #ln[6] phi0 - phase. Must be zero.
                I.Description = comment
                FF.ImproperTypes.append(copy(I))
        elif field=='HBOND':
            continue

        elif field=='END':
            continue
            #break
        else:
            continue

    ifile.close()
    
    return FF






def ReadCharmmStreamMolecularTopologyFile(STRfn):

    def ReadFFVersion(STRfn):
        ifile=open(STRfn,'r')
        #for line in ifile:
        FFName='cgenff2b7'    
        ifile.close()
        return FFName
    
    def ReadProgVersion(STRfn):
        ifile=open(STRfn,'r')
        #for line in ifile:
        Prog = 'paramchem'    
        ifile.close()
        return Prog
    
    tmp = STRfn.split('.');
    MoleculeName = tmp[0]
    Mol = Molecule(MoleculeName)
    Mol.FFName = ReadFFVersion(STRfn)    
    Mol.GeneratedBy = ReadProgVersion(STRfn)
    UnitsIn={ 'Energy':'kcal', 'Distance':'angstr', 'Angle':'rad', 'Number':'mol' }
    Mol.Units = UnitsIn

    # Due to AUTOGENERATE  ANGLE DIHEDRAL, there are usually no ANGLE and DIHEDRAL entries in the residue topology file in Charmm toppar:
    # BUT not always!

    FieldsList=['ATOM','BOND','ANGLE','PROPER','IMPROPER']
    ifile=open(STRfn,'r')

    iAt = 0
    iBd = 0
    iAn = 0
    iTn = 0
    iIm = 0


 
    while 1:
        
        (line,comment,lnextline,lexit) = GetLineCharmmTopParFile(ifile)        
        if lexit: break
        if lnextline: continue
        ln=line.split()

        #print '!!!!! ',ln,lexit,lnextline
        if ln[0][0:4]=='RESI':
            # Start reading data for the residue toology until the key END
            MoleculeName=ln[1]
        
            field=ln[0]
            while field != 'END':
                (line,comment,lnextline,lexit) = GetLineCharmmTopParFile(ifile)        
                if lexit: break
                if lnextline: continue

                ln=line.split()
                field=ln[0]
            
#ATOM O1     OG2D1  -0.524 !    0.000
#ATOM C1     CG2O1   0.527 !    0.000
                if field[:4]=='ATOM':
                    iAt += 1
                    A = Atom(iAt)
                    A.Name = ln[1] # Mind! atom names are unique within given residue
                    A.Type = ln[2]
                    A.Charge = float(ln[3])
                    A.Description = "Qpenalty = "+comment
                    A.Func = 'LJ_Rmin2_eps'
                    Mol.Atoms.append(copy(A))

#BOND H7   C8  
#BOND H9   C8   H5   N1  
#DOUBLE H O
#TRIPLE H O
                if ((field[0:4]=='BOND') or (field[0:4]=='DOUB') or (field[0:4]=='TRIP')):
                #if field[0:4]=='BOND' or field[0:4]=='DOUB':
                    for i in range(  (len(ln)-1)/2  ):
                        iBd += 1
                        B = Bond(iBd) 
                        A1=ln[2*i+1]; A2=ln[2*i+2]
                        # Find the indexes based on atom names (which are unique!)
                        
                        Atom1 = FindObjByFieldValue(Mol.Atoms,'Name',A1)
                        Atom2 = FindObjByFieldValue(Mol.Atoms,'Name',A2)
                        B.AtNumbers = [ Atom1.Number , Atom2.Number ] 
                        B.AtNames = [ A1, A2 ]
                        if field[:4]=='DOUB':
                            comment+=' : type=DOUBLE'
                        if field[:4]=='TRIP':
                            comment+=' : type=TRIPLE'
                        B.Description = comment
                        B.Func = 'Harmonic'
                        Mol.Bonds.append(copy(B))

#IMPR C1     C8     N1     O1    
                elif field[:4]=='IMPR':
                    for i in range(  (len(ln)-1)/4  ):
                        iIm += 1
                        I = Improper(iIm) 
                
                        A1=ln[4*i+1]; A2=ln[4*i+2]; A3=ln[4*i+3]; A4=ln[4*i+4]  
                        # Find the indexes based on atom names (which are unique!)
                        Atom1 = FindObjByFieldValue(Mol.Atoms,'Name',A1)
                        Atom2 = FindObjByFieldValue(Mol.Atoms,'Name',A2)
                        Atom3 = FindObjByFieldValue(Mol.Atoms,'Name',A3)
                        Atom4 = FindObjByFieldValue(Mol.Atoms,'Name',A4)
                        I.AtNumbers = [ Atom1.Number , Atom2.Number , Atom3.Number, Atom4.Number ]
                        I.AtNames = [ A1, A2, A3, A4 ]
                        I.Description = comment
                        I.Func = 'Harmonic'
                        Mol.Impropers.append(copy(I))
        
    ### Read in extra bonded FF types from stream file
    ifile.close()
    FF = ReadCharmmForceFieldFile(STRfn)
    FF.Comment += " Mind!!! These parameters were generated by analogy by paramchem server. Its validity must be checked. Please check the penalties at the FF description. ;!--- Penalties lower than 10 indicate the analogy is fair; penalties between 10   ;!--- and 50 mean some basic validation is recommended; penalties higher than   ;!--- 50 indicate poor analogy and mandate extensive validation/optimization."

    # Add missing data 
    Mol.GenAnglesTorsionsFromBondConnectivity()
    for F in Mol.Angles:
        F.Func = 'UreyBradley'
    for F in Mol.Torsions:
        F.Func = 'PeriodicMultiple'

    return Mol,FF








def GetLineCharmmTopParFile(ifile):
        
        HyphenAtTheEndOfLine=True
        lnextline=False; lexit=False
        comment=''
        cnt=-1
        while HyphenAtTheEndOfLine:
            cnt+=1
            if cnt>0:
                print 'Note!!! While reading charmm toppar file a line has " -" at the end. Skipping comments in the lines and merging the next line with the previous one.'
    
            line=ifile.readline()
            if not line: 
                lexit=True
                break
        
            # get rid of end of line character
            line=line[:-1]
        
        
            if '!' in line: 
                ln=line.split('!',1)
                line=ln[0]
                try:
                    comment=ln[1]
                except:
                    comment='no'
            
        
            if len(line)==0: 
                line=''; comment=''; lnextline=True; lexit=False
                return line,comment,lnextline,lexit
        
            if line[0]=='*': 
                try:
                    comment=line[1:-1]
                except:
                    comment=''
                line=''; comment=comment; lnextline=True; lexit=False
                return line,comment,lnextline,lexit
       
            if line[-1]=='-': 
                # Merge lines when there is "-" at the end of line.
                # Mind!!! Now only one line can be merged.
                HyphenAtTheEndOfLine=True
            else:
                HyphenAtTheEndOfLine=False
                
        return line,comment,lnextline,lexit


