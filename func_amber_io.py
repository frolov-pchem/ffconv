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
import numpy as np
NaN=float('NaN')

###
### Set global variables
### 
user=str(getpass.getuser())
#today=datetime.today()
today=str(today.strftime('%m/%d/%Y  %H:%M:%S'))




def GetArrayFromAmberTopology(lines,field):
# The fuction takes a field name, reads appropriat format of the data, and fill the data into an array. 
# Gets an array with lines of the file. field - filed name.
    OutArr=[]
    for ind in range(len(lines)):
        line=lines[ind]
        ln=line.split()
        if (len(ln)<=1): continue
        if (ln[1] == field):
            # Finding the length of entries in the array
            fmt=lines[ind+1][:-1]; fmt=fmt.strip(); fmt=fmt[:-1]; fmt=fmt.split('('); fmt=fmt[1]; fmt=fmt.split('.'); fmt=fmt[0]
            for i in range(len(fmt)):
                try:
                    float(fmt[i])
                except:
                    break
            fmt=int(fmt[i+1:])
            #print field,"format", fmt


            # Filling in the array
            tmpind=ind+2
            while lines[tmpind][0]!="%":
                line2 = lines[tmpind]; l=len(line2[:-1])
                for j in range(0,l,fmt):
                     OutArr.append(line2[j:j+fmt])
                tmpind+=1
                if tmpind > len(lines)-1:
                    break
    return OutArr


def PrintAmberTopFieldValues(lst,step,fmt):
    i=0
    s=''
    while i < len(lst):
        step = min(step,len(lst)-i)
        fmt2 = step * fmt+"\n"
        s += fmt2 % tuple(lst[i:i+step])  
        i+=step
    return s







def ReadAmberTopologyFile(fn):

### Ignore the following flags in the Amber topology
#%FLAG NUMBER_EXCLUDED_ATOMS                                                     
#%FLAG NONBONDED_PARM_INDEX                                                      
#%FLAG RESIDUE_POINTER
#%FLAG HBOND_ACOEF                                                               
#%FLAG HBOND_BCOEF                                                               
#%FLAG HBCUT                                                                     
#%FLAG JOIN_ARRAY                                                                
#%FLAG IROTAT                                                                    
#%FLAG SOLTY                                                                     
#%FLAG EXCLUDED_ATOMS_LIST                                                       

# Add these arrays to Mol.ExtraData
#%FLAG TREE_CHAIN_CLASSIFICATION                                                 
#%FLAG RADIUS_SET                                                                
#%FLAG RADII                                                                     
#%FLAG SCREEN                                                                    

    ### Constants
    FFName = 'GAFF'    
    # Implied Combining rule: Lorentz-Bershelot
    # => we convert the A and B LJ parameters to Rmin/2 and eps (as they originally present in Amber FFs)
    CombRule='LorentzBerthelot'    

    AtomList = [] ; AtType=[] ; AtName=[] ; AtFunc = [] ; q=[]
    BondList = [] ; BdTypes=[] ; BdNames = [] ; BdNumbers=[] ; BdFunc = []
    PairList = [] ; PrTypes=[] ; PrNames = [] ; PrNumbers=[] ; PrFunc = []
    AngleList = [] ; AnTypes=[] ; AnNames = [] ; AnNumbers=[] ; AnFunc = []
    TorsionList = [] ; TnTypes=[] ; TnNames = [] ; TnNumbers=[] ; TnFunc = []
    ImproperList = [] ; ImTypes = [] ; ImNames = [] ; ImNumbers=[] ; ImFunc = []
    pAt = []; pAt_14 = [] ; pBd = []; pAn = []; pTn = [] ; pIm = [] ; pPr = []

    FFAtType=[] ; FFAtFunc = [] ; FFq=[] ; FFpAt=[] ; FFpAt_14=[]
    FFBdTypes=[] ; FFBdNames = [] ; FFBdFunc = []
    FFPrTypes=[] ; FFPrNames = [] ; FFPrFunc = []
    FFAnTypes=[] ; FFAnNames = [] ; FFAnFunc = []
    FFTnTypes=[] ; FFTnNames = [] ; FFTnFunc = []
    FFImTypes = [] ; FFImNames = [] ; FFImFunc = []
    FFpBd = [] ; FFpPr = [] ; FFpAn = [] ; FFpTn = [] ; FFpIm = []

    
    ifile=open(fn,'r')
    lines=ifile.readlines()
    
    # READ IN NESESSARY DATA
    ATOM_NAME = GetArrayFromAmberTopology(lines,'ATOM_NAME')
    CHARGE = GetArrayFromAmberTopology(lines,'CHARGE')
    ATOM_TYPE_INDEX = GetArrayFromAmberTopology(lines,'ATOM_TYPE_INDEX')
    LENNARD_JONES_ACOEF = GetArrayFromAmberTopology(lines,'LENNARD_JONES_ACOEF')
    LENNARD_JONES_BCOEF = GetArrayFromAmberTopology(lines,'LENNARD_JONES_BCOEF')
    AMBER_ATOM_TYPE = GetArrayFromAmberTopology(lines,'AMBER_ATOM_TYPE')
    MASS = GetArrayFromAmberTopology(lines,'MASS')

    BOND_FORCE_CONSTANT  = GetArrayFromAmberTopology(lines,'BOND_FORCE_CONSTANT')
    BOND_EQUIL_VALUE  = GetArrayFromAmberTopology(lines,'BOND_EQUIL_VALUE')
    ANGLE_FORCE_CONSTANT  = GetArrayFromAmberTopology(lines,'ANGLE_FORCE_CONSTANT')                                                     
    ANGLE_EQUIL_VALUE    = GetArrayFromAmberTopology(lines,'ANGLE_EQUIL_VALUE')                                                      
    DIHEDRAL_FORCE_CONSTANT  = GetArrayFromAmberTopology(lines,'DIHEDRAL_FORCE_CONSTANT')                                                     
    DIHEDRAL_PERIODICITY  = GetArrayFromAmberTopology(lines,'DIHEDRAL_PERIODICITY')
    DIHEDRAL_PHASE  = GetArrayFromAmberTopology(lines,'DIHEDRAL_PHASE')
    BONDS_INC_HYDROGEN  = GetArrayFromAmberTopology(lines,'BONDS_INC_HYDROGEN')
    BONDS_WITHOUT_HYDROGEN  = GetArrayFromAmberTopology(lines,'BONDS_WITHOUT_HYDROGEN')
    ANGLES_INC_HYDROGEN  = GetArrayFromAmberTopology(lines,'ANGLES_INC_HYDROGEN')
    ANGLES_WITHOUT_HYDROGEN  = GetArrayFromAmberTopology(lines,'ANGLES_WITHOUT_HYDROGEN')
    DIHEDRALS_INC_HYDROGEN  = GetArrayFromAmberTopology(lines,'DIHEDRALS_INC_HYDROGEN')
    DIHEDRALS_WITHOUT_HYDROGEN  = GetArrayFromAmberTopology(lines,'DIHEDRALS_WITHOUT_HYDROGEN')
    
    SCEE_SCALE_FACTOR  = GetArrayFromAmberTopology(lines,'SCEE_SCALE_FACTOR')
    tmp=np.array(SCEE_SCALE_FACTOR,dtype='f'); tmp = tmp[tmp!=0]
    if len(tmp)!=0:
        fudgeQQ = CheckIfAllElementsInArrayAreTheSameAndReturn(tmp)
        fudgeQQ = 1.0/fudgeQQ    
    else:
        # It seems that there is no SCEE_SCALE_FACTOR in the Amber topology file (maybe old version fromat ? )
        # set fudgeQQ to default value for GAFF
        fudgeQQ = 0.83333333333333

    
    SCNB_SCALE_FACTOR  = GetArrayFromAmberTopology(lines,'SCNB_SCALE_FACTOR')
    tmp=np.array(SCNB_SCALE_FACTOR,dtype='f'); tmp = tmp[tmp!=0]
    if len(tmp)!=0:
        fudgeLJ = CheckIfAllElementsInArrayAreTheSameAndReturn(tmp)
        fudgeLJ = 1.0/fudgeLJ 
    else:
        # It seems that there is no SCNB_SCALE_FACTOR in the Amber topology file (maybe old version fromat ? )
        # set fudgeLJ to default value for GAFF
        fudgeLJ = 0.5
        

    RESIDUE_LABEL  = GetArrayFromAmberTopology(lines,'RESIDUE_LABEL')
    ResName = CheckIfAllElementsInArrayAreTheSameAndReturn(RESIDUE_LABEL)

    tmp=fn.split('/'); tmp=tmp[-1]; tmp=tmp.split('.'); name=tmp[0]
    MoleculeName = name



    #### Atoms
    NAtom=len(ATOM_TYPE_INDEX)
    for val in ATOM_NAME:
        AtName.append(val.strip())

    for i,val in enumerate(ATOM_TYPE_INDEX):
        ATOM_TYPE_INDEX[i] = val.strip()

    for i,val in enumerate(AMBER_ATOM_TYPE):
        AMBER_ATOM_TYPE[i] = val.strip()

    NBAtTypeIndex = ATOM_TYPE_INDEX
    AtTypeNames = AMBER_ATOM_TYPE
    AtType = AtTypeNames
    
    for i,val in enumerate(MASS):
        MASS[i] = float(val)
    AtMass = MASS
    

    #lst=[]
    map_AT_NBAT={}
    lst=[]
    for i,v in enumerate(AtTypeNames):
        if not v in lst:
            lst.append(v)
            NBind = int(float(NBAtTypeIndex[i]))
            #if not v in map_AT_NBAT.keys()
            #    map_AT_NBAT[NBind] = [ v ]
            #else:
            #    map_AT_NBAT[NBind].append(v)
            map_AT_NBAT[v] = NBind   
            

 
    AtomList = [i+1 for i in range(NAtom)]
        
    for i in range(NAtom):
        ind=int(ATOM_TYPE_INDEX[i])
        q.append(float(CHARGE[i])/18.2223)
    
    Aij = LENNARD_JONES_ACOEF
    Bij = LENNARD_JONES_BCOEF
    rmin2ij=[]; epsij=[]
    for j in range(len(Aij)):
        if (float(Bij[j]) != 0 and float(Aij[j] != 0)):
            rmin2ij.append( 2**0.1666666666666666/2*(float(Aij[j])/float(Bij[j]))**0.166666666666666666666)
            epsij.append((float(Bij[j])**2/4/float(Aij[j])))
        else:
            rmin2ij.append(0.0)
            epsij.append(0.0)

    # CHOOSING ONLY sigii and epsii terms, throwing away ij terms. ONLY ATOM TYPES!!!
    NNBAtTy=len(list(set(ATOM_TYPE_INDEX)))
    cnt=-1;
    ki=0
    map_NBATind_ATind={}

    D={}
    tmp_FFpAt = []
    tmp_FFpAt_14 = []
    for i in range( NNBAtTy ):
        for j in range( i+1 ):
            cnt+=1
        
        #for k,v in map_AT_NBAT.iteritems():
            #print k,v
            #if v == i+1: 
            #     ki += 1
        tmp_FFpAt.append( [ rmin2ij[cnt], epsij[cnt] ] )
        tmp_FFpAt_14.append( [ rmin2ij[cnt], epsij[cnt]*fudgeLJ ] )
                #map_NBATind_ATind[i+1] = ki

    #print map_AT_NBAT
    #print map_NBATind_ATind
        #FFAtFunc.append( copy(AtFunc[cnt]) )
    
    
    # ASSIGNING PARAMETERS FROM ATOMTYPES TO ATOMS
    #NAtom=len(ATOM_TYPE_INDEX)
    #AtType=[]
    FFAtMass = []
    for i in range(NAtom):

        AtTy = AtTypeNames[i]         
        NBATind = map_AT_NBAT[AtTy] 
        #ind = map_NBATind_ATind[NBATind]
        ind = NBATind
        pAt.append( copy(tmp_FFpAt[ind-1]) )
        pAt_14.append( copy(tmp_FFpAt_14[ind-1]) )
        AtFunc.append( 'LJ_Rmin2_eps' )

        #print 'type, NBATind, realATind  ',AtTy, NBATind, ind
        #AtType.append(AtTy)

        if not AtType[i] in FFAtType:
            FFAtType.append( copy(AtType[i]) )
            FFAtFunc.append( copy(AtFunc[i]) )
            FFpAt.append( copy(pAt[i]) )
            FFpAt_14.append( copy(pAt[i]) )
            FFAtMass.append( copy(AtMass[i]) )
   
 
    #print list(enumerate(FFpAt))
    #print list(enumerate(FFAtType))
    #stop


    #### Bonds
    pBdTypes = np.c_[  np.array(BOND_EQUIL_VALUE,dtype='f'), np.array(BOND_FORCE_CONSTANT,dtype='f') * 2.0   ] # Mind conversion of force constant from K to K/2 between Amber and internal convention
    
    arr = np.array( BONDS_WITHOUT_HYDROGEN + BONDS_INC_HYDROGEN,dtype='i' )
    arr = np.reshape(arr, (len(arr)/3,3) )
    arr = np.c_[ arr[:,0]/3+1, arr[:,1]/3+1, arr[:,2] ] # change the pointers for coordinate array to real atom indexes
    for i in range(len(arr)):
        BdTypes.append( [ AtType[ arr[i,0]-1 ], AtType[ arr[i,1]-1 ] ]  ) 
        BdNames.append( [ AtName[ arr[i,0]-1 ], AtName[ arr[i,1]-1 ] ]  ) 
        BdNumbers.append( [ arr[i,0], arr[i,1] ] ) 
        #print arr[i,2]-1
        pBd.append( list(pBdTypes[ arr[i,2]-1 ,: ]) )
        BdFunc.append( 'Harmonic' )
        
        tp=[ AtType[ arr[i,0]-1 ], AtType[ arr[i,1]-1 ] ]
        if not tp in FFBdTypes and not tp[::-1] in FFBdTypes:
            # Add BondType
            FFBdTypes.append( copy( BdTypes[-1] ) )
            FFpBd.append( copy( pBd[-1] ) )
            FFBdFunc.append( copy( BdFunc[-1] ) )
            

    #### Angles
    pAnTypes = np.c_[  np.array(ANGLE_EQUIL_VALUE,dtype='f'), np.array(ANGLE_FORCE_CONSTANT,dtype='f') * 2.0  ] # Mind conversion of force constant from K to K/2 between Amber and internal convention
    
    arr = np.array( ANGLES_WITHOUT_HYDROGEN + ANGLES_INC_HYDROGEN , dtype='i' )
    arr = np.reshape(arr, (len(arr)/4,4) )
    arr = np.c_[ arr[:,0]/3+1, arr[:,1]/3+1, arr[:,2]/3+1 , arr[:,3] ] # change the pointers for coordinate array to real atom indexes
    for i in range(len(arr)):
        AnTypes.append( [ AtType[ arr[i,0]-1 ], AtType[ arr[i,1]-1 ], AtType[ arr[i,2]-1 ] ] ) 
        AnNames.append( [ AtName[ arr[i,0]-1 ], AtName[ arr[i,1]-1 ], AtName[ arr[i,2]-1 ] ]  ) 
        AnNumbers.append( [ arr[i,0], arr[i,1], arr[i,2] ] ) 
        pAn.append( list(pAnTypes[ arr[i,3]-1 ,:]) )
        AnFunc.append( 'Harmonic' )
        
        tp=[ AtType[ arr[i,0]-1 ], AtType[arr[i,1]-1], AtType[ arr[i,2]-1 ] ]
        if not tp in FFAnTypes and not tp[::-1] in FFAnTypes:
            FFAnTypes.append( copy( AnTypes[-1] ) )
            FFpAn.append( copy( pAn[-1] ) )
            FFAnFunc.append( copy( AnFunc[-1] ) )
            
    
    #### Torsions, Impropers, Pairs14
    #pTnTypes = np.c_[  np.array(DIHEDRAL_PHASE,dtype='f'), np.array(DIHEDRAL_FORCE_CONSTANT,dtype='f') / 2.0, np.array(DIHEDRAL_PERIODICITY,dtype='f') ] # Mind conversion to Amber notation.

    ### TODO: 
    # it seems that potential parameters [Vn] for dihedrals in amber topology do not match the functional form [ Vn/2 * (1 + cos( n*phi - gamma ))].
    # the given potential parameters are for this func form: [ Vn * (1 + cos( n*phi - gamma ))]
    pTnTypes = np.c_[  np.array(DIHEDRAL_PHASE,dtype='f'), np.array(DIHEDRAL_FORCE_CONSTANT,dtype='f'), np.array(DIHEDRAL_PERIODICITY,dtype='f') ] # Mind conversion to Amber notation.
    dih_list = DIHEDRALS_WITHOUT_HYDROGEN + DIHEDRALS_INC_HYDROGEN
    if '-0' in dih_list:
        # Reverting the dihedral where -0
        i=dih_list.index('-0')/5
        dih_list[i:i+4] = reversed(dih_list[i:i+4])  
        if '-' in dih_list[i+2]: 
            dih_list[i+1].replace('-','')
            dih_list[i+2]='-'+dih_list[i+2]
        if '-' in dih_list[i+3]: 
            dih_list[i].replace('-','')
            dih_list[i+3]='-'+dih_list[i+3]
    
    
    arr = np.array( dih_list, dtype='i' )
    arr = np.reshape(arr, (len(arr)/5,5) )
    ndx=np.zeros( (int(np.size(arr,0)) ,int(np.size(arr,1))-1), dtype='i' ) 
    for i in range(len(ndx)):
        for j in range(np.size(ndx,1)):
            ndx[i,j] = int(math.modf( np.abs(arr[i,j])/3)[1])+1
    
    #arr = np.c_[ (arr[:,0]/3)[1])+1,int(math.modf(arr[:,1]/3)[1])+1,int(math.modf(arr[:,2]/3)[1])+1,int(math.modf(arr[:,3]/3)[1])+1, arr[:,4] ] # change the pointers for coordinate array to real atom indexes
    for i in range(len(arr)):
        tp=[ AtType[ ndx[i,0]-1 ], AtType[ ndx[i,1]-1 ], AtType[ ndx[i,2]-1 ], AtType[ ndx[i,3]-1 ] ]
        tn=[ AtName[ ndx[i,0]-1 ], AtName[ ndx[i,1]-1 ], AtName[ ndx[i,2]-1 ], AtName[ ndx[i,3]-1 ] ]
        if arr[i,3] >= 0:
            TnTypes.append( tp ) 
            TnNames.append( tn ) 
            TnNumbers.append( [ ndx[i,0], ndx[i,1], ndx[i,2], ndx[i,3] ] ) 
            pTn.append( list(pTnTypes[ arr[i,4]-1 ,:])  )
            TnFunc.append( 'PeriodicMultiple' )
            #print FFTnTypes
            #print FFpTn
            
            l11=(tp in FFTnTypes)
            l21=(tp[::-1] in FFTnTypes)
            
            l12=False
            for ii in range(len(FFTnTypes)):
                if FFTnTypes[ii] == tp:
                    l12=(pTn[-1] == FFpTn[ii])
                    if l12: break
            
            l22=False
            for ii in range(len(FFTnTypes)):
                if FFTnTypes[ii] == tp[::-1]:
                    l22=(pTn[-1] == FFpTn[ii])
                    if l22: break
            
            if not (l11 and l12) and not (l21 and l22):
                FFTnTypes.append( copy( TnTypes[-1] ) )
                FFpTn.append( copy( pTn[-1] ) )
                FFTnFunc.append( copy( TnFunc[-1] ) )
        else:
            ImTypes.append( tp ) 
            ImNames.append( tn ) 
            ImNumbers.append( [ ndx[i,0], ndx[i,1], ndx[i,2], ndx[i,3] ] ) 
            pIm.append( list(pTnTypes[ arr[i,4]-1 ,:])  )
            ImFunc.append( 'Periodic' )

            l11=(tp in FFImTypes)
            l21=(tp[::-1] in FFImTypes)
            
            l12=False
            for ii in range(len(FFImTypes)):
                if FFImTypes[ii] == tp:
                    l12=(pIm[-1] == FFpIm[ii])
                    if l12: break
            
            l22=False
            for ii in range(len(FFImTypes)):
                if FFImTypes[ii] == tp[::-1]:
                    l22=(pIm[-1] == FFpIm[ii])
                    if l22: break
            
            if not (l11 and l12) and not (l21 and l22):
                FFImTypes.append( copy( ImTypes[-1] ) )
                FFpIm.append( copy( pIm[-1] ) )
                FFImFunc.append( copy( ImFunc[-1] ) )


            
        if arr[i,2] >= 0:
            # Pair 1-4 should be taken into account
            PrTypes.append( [ AtType[ ndx[i,0]-1 ], AtType[ ndx[i,3]-1 ] ]  ) 
            PrNames.append( [ AtName[ ndx[i,0]-1 ], AtName[ ndx[i,3]-1 ] ]  ) 
            PrNumbers.append( [ ndx[i,0], ndx[i,3] ] ) 
            pPr.append( ApplyCombRules( pAt_14[ndx[i,0]-1], pAt_14[ndx[i,3]-1], CombRule, 'LJ_Rmin2_eps' ) )
            PrFunc.append( 'LJ_Rmin2_eps' )
            if not [ AtType[ ndx[i,0]-1 ], AtType[ ndx[i,3]-1 ] ] in FFPrTypes and not [ AtType[ ndx[i,3]-1 ], AtType[ ndx[i,0]-1 ] ] in FFPrTypes:
                FFPrTypes.append( copy( PrTypes[-1] ) )
                FFpPr.append( copy( pPr[-1] ) )
                FFPrFunc.append( copy( PrFunc[-1] ) )
        else:
            pass
    

    ifile.close()
   


    UnitsIn={ 'Energy':'kcal', 'Distance':'angstr', 'Angle':'rad', 'Number':'mol' }

    ###
    ### Fill in Molecule
    ###
    Mol = Molecule(MoleculeName)
    Mol.FFName = FFName
    Mol.Units = UnitsIn    
    Mol.GeneratedBy='antechamber+tleap'

    #FillTopParAttributes(Mol, 'Atom', Name=AtName, Number=AtomList, Type=AtType, Charge=q, Func=AtFunc, p=pAt, p_14=pAt_14, Mass=[], PeriodicTableNum=[], X=[], Y=[], Z=[])
    for i,v in enumerate(AtName):
        F = Atom(i+1)
        F.Name = AtName[i]
        F.Type = AtType[i]
        F.Charge = q[i]
        F.Func = AtFunc[i]
        F.p = pAt[i]
        F.p_14 = pAt_14[i]
        F.Mass = AtMass[i]
        Mol.Atoms.append(copy(F))
        
    #FillTopParAttributes(Mol, 'Bond', AtTypes=BdTypes, AtNames=BdNames, AtNumbers=BdNumbers, Func=BdFunc, p=pBd)
    for i,v in enumerate(BdTypes):
        F = Bond(i+1)
        F.AtNames = BdNames[i]
        F.AtTypes = BdTypes[i]
        F.AtNumbers = BdNumbers[i]
        F.Func = BdFunc[i]
        F.p = pBd[i]
        Mol.Bonds.append(copy(F))
        
    #FillTopParAttributes(Mol, 'Pair', AtTypes=PrTypes, AtNames=PrNames, AtNumbers=PrNumbers, Func=PrFunc, p=pPr)
    for i,v in enumerate(PrTypes):
        F = Pair(i+1)
        F.AtNames = PrNames[i]
        F.AtTypes = PrTypes[i]
        F.AtNumbers = PrNumbers[i]
        F.Func = PrFunc[i]
        F.p = pPr[i]
        Mol.Pairs.append(copy(F))

    #FillTopParAttributes(Mol, 'Angle', AtTypes=AnTypes, AtNames=AnNames, AtNumbers=AnNumbers, Func=AnFunc, p=pAn)
    for i,v in enumerate(AnTypes):
        F = Angle(i+1)
        F.AtNames = AnNames[i]
        F.AtTypes = AnTypes[i]
        F.AtNumbers = AnNumbers[i]
        F.Func = AnFunc[i]
        F.p = pAn[i]
        Mol.Angles.append(copy(F))

    #FillTopParAttributes(Mol, 'Torsion', AtTypes=TnTypes, AtNames=TnNames, AtNumbers=TnNumbers, Func=TnFunc, p=pTn)
    for i,v in enumerate(TnTypes):
        F = Torsion(i+1)
        F.AtNames = TnNames[i]
        F.AtTypes = TnTypes[i]
        F.AtNumbers = TnNumbers[i]
        F.Func = TnFunc[i]
        F.p = pTn[i]
        Mol.Torsions.append(copy(F))

    #FillTopParAttributes(Mol, 'Improper', AtTypes=ImTypes, AtNames=ImNames, AtNumbers=ImNumbers, Func=ImFunc, p=pIm)
    for i,v in enumerate(ImTypes):
        F = Improper(i+1)
        F.AtNames = ImNames[i]
        F.AtTypes = ImTypes[i]
        F.AtNumbers = ImNumbers[i]
        F.Func = ImFunc[i]
        F.p = pIm[i]
        Mol.Impropers.append(copy(F))

#%FLAG TREE_CHAIN_CLASSIFICATION                                                 
#%FLAG RADIUS_SET                                                                
#%FLAG RADII                                                                     
#%FLAG SCREEN                                                                    
    TREE_CHAIN_CLASSIFICATION = GetArrayFromAmberTopology(lines,'TREE_CHAIN_CLASSIFICATION')
    RADIUS_SET = GetArrayFromAmberTopology(lines,'RADIUS_SET')
    RADII = GetArrayFromAmberTopology(lines,'RADII')
    SCREEN = GetArrayFromAmberTopology(lines,'SCREEN')
    Mol.ExtraData['TREE_CHAIN_CLASSIFICATION']=TREE_CHAIN_CLASSIFICATION 
    Mol.ExtraData['BORN_RADIUS_SET']=RADIUS_SET
    Mol.ExtraData['BORN_RADII']=np.array(RADII,dtype='f')
    Mol.ExtraData['BORN_SCREEN']=np.array(SCREEN,dtype='f')
    
    #np.array(BOND_EQUIL_VALUE,dtype='f')


    ###
    ### Fill in FF
    ###
    # TODO get FFname from cmdline arg
    FF = ForceField('gaff')
    FF.Units = UnitsIn
    FF.NREXCL = 3
    FF.GeneratedBy = 'antechamber+tleap'
    FF.Comment = ''

    #FillTopParAttributes(FF, 'AtomType', Type=FFAtType, Charge=FFq, Func=FFAtFunc, p=FFpAt, p_14=FFpAt_14, Mass=[], PeriodicTableNum=[])
    for i,v in enumerate(FFAtType):
        F = AtomType()
        F.Type = FFAtType[i]
        F.Charge = 0.0 # dummy atom charge
        F.Func = FFAtFunc[i]
        F.p = FFpAt[i]
        F.p_14 = FFpAt_14[i]
        F.Mass = FFAtMass[i]
        FF.AtomTypes.append(copy(F))

    #FillTopParAttributes(FF, 'BondType', AtTypes=FFBdTypes, Func=FFBdFunc, p=FFpBd)
    for i,v in enumerate(FFBdTypes):
        F = BondType()
        F.AtTypes = FFBdTypes[i]
        F.Func = FFBdFunc[i]
        F.p = FFpBd[i]
        FF.BondTypes.append(copy(F))

    #FillTopParAttributes(FF, 'PairType', AtTypes=FFPrTypes, Func=FFPrFunc, p=FFpPr)
    for i,v in enumerate(FFPrTypes):
        F = PairType()
        F.AtTypes = FFPrTypes[i]
        F.Func = FFPrFunc[i]
        F.p = FFpPr[i]
        FF.PairTypes.append(copy(F))

    #FillTopParAttributes(FF, 'AngleType', AtTypes=FFAnTypes, Func=FFAnFunc, p=FFpAn)
    for i,v in enumerate(FFAnTypes):
        F = AngleType()
        F.AtTypes = FFAnTypes[i]
        F.Func = FFAnFunc[i]
        F.p = FFpAn[i]
        FF.AngleTypes.append(copy(F))

    #FillTopParAttributes(FF, 'TorsionType', AtTypes=FFTnTypes, Func=FFTnFunc, p=FFpTn)
    for i,v in enumerate(FFTnTypes):
        F = TorsionType()
        F.AtTypes = FFTnTypes[i]
        F.Func = FFTnFunc[i]
        F.p = FFpTn[i]
        FF.TorsionTypes.append(copy(F))

    #FillTopParAttributes(FF, 'ImproperType', AtTypes=FFImTypes, Func=FFImFunc, p=FFpIm)
    for i,v in enumerate(FFImTypes):
        F = ImproperType()
        F.AtTypes = FFImTypes[i]
        F.Func = FFImFunc[i]
        F.p = FFpIm[i]
        FF.ImproperTypes.append(copy(F))

    
    #print Mol
    #print FF
    #stop
    return Mol,FF





#----------------------------------------------------------------------------------------------------------------------
def WriteAmberMoleculeTopFile(FName,Mol,FF,**kwargs):
#
#   EXAMPLE USAGE
#
#   DESCRIPTION
#   
#   INPUT
# FName - name of the output file (string)
# Mol - a Molecule class object
# FF - a ForceField class object
# (Optional)
#
#   REQUIREMENTS
#
#
    ### TODO: 
   

    ###
    ### Important info: NAB does not support:
    ### 1) Chamber type of Amber topology
    ### 2) SCEE and SCLJ 
    ### AND DOES NOT WRITE ANY ERROR MESSAGES !
    ###
 
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
    ### Defaults
    ### 
    Vars={}
    
    # Set default values for Vars:
    Vars['lBonded']='no'
    Vars['VERSION_STAMP']="V0001.000"
    Vars['DATE'] = today
    Vars['lChamber'] = 'no' # CHAMBER
    Vars['NREXCL'] = FF.NREXCL # CHAMBER


    FFFamily = GetFFFamily(Mol.FFName)
    #if FFFamily == 'CHARMM':
    #    Vars['lChamber'] = 'yes'
    #else:
    #    Vars['lChamber'] = 'no'

    ###
    ### Process function arguments 
    ###
    for k,v in kwargs.iteritems():
        if not k in Vars.keys():
            sys.stderr.write('--- !Warning in function '+inspect.stack()[0][3]+' not known key argument '+k+'. Continue.')
        for key,val in Vars.iteritems():
            if k==key:
                Vars[k]=v
            else:
                pass

    lBonded = Vars['lBonded']
    VERSION_STAMP = Vars['VERSION_STAMP']
    DATE=Vars['DATE']
    lChamber=Vars['lChamber']
    NREXCL=Vars['NREXCL']
    
    #NREXCL=0    

    ###
    ### Checks
    ###

        
    
    ###
    ### Convert NB params: e.g from LJ_rmin2_eps tp LJ_sig_eps
    ###
    #ConvertNonBondedParamsTo(Mol,NBFuncType)

    ###
    ### Perform conversions
    ###

    # Generate exclusions
    Mol.GenExclusions(NREXCL)
    
    # Convert units to AMBER
    ConvertUnitsTo(Mol,'toamber')

    # Convert angle units to degree 
    #ConvertUnitsTo(Mol,'toAMBER',lrad2degree='yes')



    ###
    ### Get missing data
    ###


    # Atom types
    ATOM_TYPE_INDEX=[]
    lst_ATN=[]
    lst_AT=[]
    cnt=0
    D_AT_ATNUM={}
    D_ATNUM_AT={}
    for A in Mol.Atoms: 
        AT = FindObjByFieldValue(FF.AtomTypes, 'Type', A.Type)        
        if not AT.Type in lst_AT:
            cnt+=1
            lst_ATN.append(  cnt  )
            lst_AT.append(AT.Type)
            D_AT_ATNUM[AT.Type] = cnt
            D_ATNUM_AT[cnt] = AT.Type
            ATOM_TYPE_INDEX.append(cnt)
        else:
            ATOM_TYPE_INDEX.append(D_AT_ATNUM[AT.Type])

    NTYPES = len(D_AT_ATNUM.keys())


    # Dihedrals
    if len(Mol.Torsions) == 0 and lBonded != "yes":
        # Maybe the 
        sys.stderr.write('--- !Note in function '+inspect.stack()[0][3]+': you requested NOT to write bonded terms to Amber topology (lBonded = ['+str(lBonded)+']). But, unfortunately, to make it possible to calculate the 1-4 interactions one must include the torsions dihedrals into topology. Now there are no Torsions in the Molecule instance: maybe they are missing? Continue.\n ')
        #Mol.GenAnglesTorsionsFromBondConnectivity()
        #sys.stderr.write('--- !Note in function '+inspect.stack()[0][3]+': ['+str(len(Mol.Torsions))+'] torsions generated. Continue.\n ')

    DIHEDRALS = copy(Mol.Torsions + Mol.Impropers)
    DIHEDRALS_H = []
    DIHEDRALS_A = []
    for D in DIHEDRALS:
        F = copy(D)
        if F.AtTypes[0][0].upper()=='H' or F.AtTypes[3][0].upper()=='H' or F.AtTypes[1][0].upper()=='H' or F.AtTypes[2][0].upper()=='H':
            DIHEDRALS_H.append(copy(F))
        else:
            DIHEDRALS_A.append(copy(F))
    cnt=0
    for F in DIHEDRALS_A: 
            for FT in FF.TorsionTypes + FF.ImproperTypes:
                if (F.AtTypes == FT.AtTypes or F.AtTypes[::-1] == FT.AtTypes) and FT.__class__.__name__[:-4] == F.__class__.__name__:
                #if F.AtTypes == FT.AtTypes or F.AtTypes[::-1] == FT.AtTypes:
                    cnt+=1
    N_DIH_TYPES_A = cnt

    cnt=0
    for F in DIHEDRALS_H: 
            for FT in FF.TorsionTypes + FF.ImproperTypes:
                if (F.AtTypes == FT.AtTypes or F.AtTypes[::-1] == FT.AtTypes) and FT.__class__.__name__[:-4] == F.__class__.__name__:
                #if F.AtTypes == FT.AtTypes or F.AtTypes[::-1] == FT.AtTypes:
                    cnt+=1
    N_DIH_TYPES_H = cnt

    N_DIH_TYPES = N_DIH_TYPES_H + N_DIH_TYPES_A


    # Exclusions
    Dexcl={}
    for F in Mol.Atoms:
        Dexcl[F.Number]=[]
    for F in Mol.Exclusions:
        Dexcl[F.AtNumbers[0]].append(F.AtNumbers[1])
    NNB=0
    for k,v in Dexcl.iteritems():
        val = len(v)
        if val == 0: val = 1   # if no exclusions for given atom set the value to 1 and then in EXCLUDED_ATOMS_LIST put zero
        NNB += val        

    ###
    ### Checks
    ###
    if lBonded == "yes":
        # TODO
        # Yet not supported. Must to extend the code
        sys.stderr.write('--- !Error in function '+inspect.stack()[0][3]+': you requested to write bonded terms to Amber topology (lBonded = ['+str(lBonded)+']), but this is not supported yet. Extend the code. Now you can write only atoms to Amber topology - which can be used e.g. in #D RISM single point calculations with NAB. The molecule instance named ['+Mol.Name+']. Aborting.\n ')
        sys.exit(1)
    else:
        sys.stderr.write('--- !Note in function '+inspect.stack()[0][3]+': requested NOT to write bonded terms to Amber topology (lBonded = ['+str(lBonded)+']) (the only supported option so far). But, unfortunately, to make it possible to calculate the 1-4 interactions one must include the torsions and improper dihedrals into topology. Continue.\n ')

        


    ###
    ### Writing the TOP file
    ###
    of=open(FName,'w')


#%VERSION  VERSION_STAMP = V0001.000  DATE = 03/19/12  12:32:35                  
    # Version
    string= "%VERSION  VERSION_STAMP = "+VERSION_STAMP+"  DATE = "+DATE+"\n"

#%FLAG TITLE                                                                     
#%FORMAT(20a4)                                                                   
#MOL                                                                             
    # Header
    if lChamber != 'yes':
        string += "%FLAG TITLE\n"
    else:
        string += "%FLAG CTITLE\n"
        
    string += "%FORMAT(20a4)\n"

#REMOVED 2019-01-30
#    string += "The top for Mol ["+Mol.Name+"] created by the ["+sys.argv[0]+"] tool, generated by user "+user+"\n"
#END REMOVED 2019-01-30

#ADDED 2019-01-30
    string += Mol.Name+"\n"
#END ADDED 2019-01-30
    
#REMOVED 2019-01-30
    # Writing FFCONV info
#    string += \
#    "-- ffconv.py info: \n" 
#    for k,v in Vars.iteritems():
#        string += ""+str(k)+" = "+str(v)+"\n"
    
    # Writing Mol info
#    string += \
#    "-- molecule parameters: \n" +\
#    "FFName = "+str(Mol.FFName)+"\n" +\
#    "Units = "+str(Mol.Units)+"\n" +\
#    "GeneratedBy = "+str(Mol.GeneratedBy)+"\n" +\
#    "Comment = "+str(Mol.Comment)+"\n" 
    
    # Writing FF info
#    string+= \
#    "-- force field parameters: \n" +\
#    "Name = "+str(FF.Name)+"\n" +\
#    "Family = "+str(FF.Family)+"\n" +\
#    "Units = "+str(FF.Units)+"\n" +\
#    "GeneratedBy = "+str(FF.GeneratedBy)+"\n" +\
#    "Comment = "+str(FF.Comment)+"\n" +\
#    "NREXCL = "+str(FF.NREXCL)+"\n" +\
#    "fudgeLJ = "+str(FF.fudgeLJ)+"\n" +\
#    "fudgeQQ = "+str(FF.fudgeQQ)+"\n" +\
#    "CombRule = "+str(FF.CombRule)+"\n" 

#END REMOVED 2019-01-30

#%FLAG POINTERS                                                                  
#%FORMAT(10I8)                                                                   
#      20       9       9      11      17      14      35      17       0       0
#      89       1      11      14      17      10      13      10      10       0
#       0       0       0       0       0       0       0       0      20       0
#       0

### Set NTYPES always equal to NATYP for simplicity!
    # NATOM    : total number of atoms 20
    NATOM = len(Mol.Atoms) 
    # NTYPES   : total number of distinct atom types 9
    NTYPES = NTYPES
    if lBonded != "yes":
        # NBONH    : number of bonds containing hydrogen 9
        NBONH = 0
        # MBONA    : number of bonds not containing hydrogen 11
        MBONA = 0 
        # NTHETH   : number of angles containing hydrogen 17
        NTHETH = 0
        # MTHETA   : number of angles not containing hydrogen 14
        MTHETA = 0

        
    # NPHIH    : number of dihedrals containing hydrogen 35
    NPHIH = N_DIH_TYPES_H
    # MPHIA    : number of dihedrals not containing hydrogen 17
    MPHIA = N_DIH_TYPES_A

    # NHPARM   : currently not used 0
    NHPARM = 0 
    # NPARM    : used to determine if addles created prmtop 0
    NPARM = 0

    # NNB      : number of excluded atom [pairs] 89
    NNB = NNB   # len(Mol.Exclusions) 
    # NRES     : number of residues 1
    NRES = 1

    if lBonded != "yes":
        # NBONA    : MBONA + number of constraint bonds 11
        NBONA = 0
        # NTHETA   : MTHETA + number of constraint angles 14
        NTHETA = 0
    
    # NPHIA    : MPHIA + number of constraint dihedrals 17
    NPHIA = N_DIH_TYPES_A


    if lBonded != "yes":
        # NUMBND   : number of unique bond types 10
        NUMBND = 0
        # NUMANG   : number of unique angle types 13
        NUMANG = 0

    # NPTRA    : number of unique dihedral types 10
    NPTRA = N_DIH_TYPES


    # NATYP    : number of atom types in parameter file, see SOLTY below 10
    NATYP = NTYPES
    # NPHB     : number of distinct 10-12 hydrogen bond pair types 0
    NPHB = 0 

    # IFPERT   : set to 1 if perturbation info is to be read in 0
    IFPERT = 0
    # NBPER    : number of bonds to be perturbed 0
    NBPER = 0
    # NGPER    : number of angles to be perturbed 0
    NGPER = 0
    # NDPER    : number of dihedrals to be perturbed 0
    NDPER = 0
    # MBPER    : number of bonds with atoms completely in perturbed group 0
    MBPER = 0 
    # MGPER    : number of angles with atoms completely in perturbed group 0
    MGPER = 0
    # MDPER    : number of dihedrals with atoms completely in perturbed groups 0
    MDPER = 0
    # IFBOX    : set to 1 if standard periodic box, 2 when truncated octahedral 0
    IFBOX = 0
    # NMXRS    : number of atoms in the largest residue 20
    NMXRS = len(Mol.Atoms)
    # IFCAP    : set to 1 if the CAP option from edit was specified 0
    IFCAP = 0
    # NUMEXTRA : number of extra points found in topology 0
    NUMEXTRA = 0
    # NCOPY    : number of PIMD slices / number of beads (TODO NOTE! Not given in the topology generated by AmberTools)
    

    string += "%FLAG POINTERS\n" 
    string += "%FORMAT(10I8)\n"  
    string += "%8i%8i%8i%8i%8i%8i%8i%8i%8i%8i\n" % (NATOM,NTYPES,NBONH,MBONA,NTHETH,MTHETA,NPHIH,MPHIA,NHPARM,NPARM)
    string += "%8i%8i%8i%8i%8i%8i%8i%8i%8i%8i\n" % (NNB,NRES,NBONA,NTHETA,NPHIA,NUMBND,NUMANG,NPTRA,NATYP,NPHB) 
    string += "%8i%8i%8i%8i%8i%8i%8i%8i%8i%8i\n" % (IFPERT,NBPER,NGPER,NDPER,MBPER,MGPER,MDPER,IFBOX,NMXRS,IFCAP)
    string += "%8i\n" % (NUMEXTRA)
     


#%FLAG ATOM_NAME                                                                 
#%FORMAT(20a4)                                                                   
#O1  C1  C8  H7  H8  H9  N1  H5  C2  C3  H1  C4  H2  C7  O2  H6  C6  H3  C5  H4  

    string += "%FLAG ATOM_NAME\n"
    string += "%FORMAT(20a4)\n"
    lst=[]
    for F in Mol.Atoms:
        lst.append(F.Name.ljust(4))
    
    s = PrintAmberTopFieldValues(lst,20,"%4s")
    
    string += s

#%FLAG CHARGE                                                                    
#%FORMAT(5E16.8)                                                                 
# -1.07652246E+01  1.21241526E+01 -3.25639790E+00  1.11757366E+00  1.11988789E+00
#  1.49909395E+00 -8.62153502E+00  5.70230434E+00  6.86251818E-01 -2.72020672E+00
#  2.34158377E+00 -3.38557578E+00  2.42823081E+00  1.75504438E+00 -9.05389553E+00
#  7.58433993E+00 -2.38380484E+00  2.80340974E+00 -2.21847391E+00  3.24324140E+00
    string+="%FLAG CHARGE\n"
    if lChamber != 'yes':
        string+="%FORMAT(5E16.8)\n"                                                                 
    else:
        string+="%FORMAT(3E24.16)\n"                                                                 
        
    lst=[]
    lst2=[]
    for F in Mol.Atoms: 
        lst.append(F.Charge * 18.2223)
        lst2.append(F.Name)

    if lChamber != 'yes':
        s = PrintAmberTopFieldValues(lst,5,"%16.8e")
    else:
        s = PrintAmberTopFieldValues(lst,3,"%24.16e") #3E24.16

        
    
    string += s



#%FLAG MASS                                                                      
#%FORMAT(5E16.8)                                                                 
#  1.60000000E+01  1.20100000E+01  1.20100000E+01  1.00800000E+00  1.00800000E+00
#  1.00800000E+00  1.40100000E+01  1.00800000E+00  1.20100000E+01  1.20100000E+01
#  1.00800000E+00  1.20100000E+01  1.00800000E+00  1.20100000E+01  1.60000000E+01
#  1.00800000E+00  1.20100000E+01  1.00800000E+00  1.20100000E+01  1.00800000E+00
    string += "%FLAG MASS\n"
    string += "%FORMAT(5E16.8)\n"
    lst=[]
    for F in Mol.Atoms: lst.append(F.Mass)
    s = PrintAmberTopFieldValues(lst,5,"%16.8e")
    string += s


#%FLAG ATOM_TYPE_INDEX                                                           
#%FORMAT(10I8)                                                                   
#       1       2       3       4       4       4       5       6       2       2
#       7       2       7       2       8       9       2       7       2       7
    string+="%FLAG ATOM_TYPE_INDEX\n"                                                           
    string+="%FORMAT(10I8)\n"                                                                   
    s = PrintAmberTopFieldValues(ATOM_TYPE_INDEX,10,"%8i")
    string += s


#%FLAG NUMBER_EXCLUDED_ATOMS                                                     
#%FORMAT(10I8)                                                                   
#       8       9       6       3       2       1       8       3       9       8
#       4       7       3       6       4       1       3       2       1       1

### My result: (seems that in Amber one exclusion is added by a mistake ? ... )
#       8       9       6       3       2       1       8       3       9       8
#       4       7       3       6       4       1       3       2       1       0
    string+="%FLAG NUMBER_EXCLUDED_ATOMS\n"                                                     
    string+="%FORMAT(10I8)\n"                                                                   
    #for k,v in Dexcl.iteritems(): 
    #    lst.append(len(v))
    #    print k,v

    Dexcl={}
    for F in Mol.Atoms:
        Dexcl[F.Number]=[]

    for F in Mol.Exclusions:
        Dexcl[F.AtNumbers[0]].append(F.AtNumbers[1])
        #Dexcl[F.AtNumbers[1]].append(F.AtNumbers[0])
    lst=[]
    for k,v in Dexcl.iteritems():
        val = len(v)
        if val == 0: val = 1   # if no exclusions for given atom set the value to 1 and then in EXCLUDED_ATOMS_LIST put zero
        lst.append(val)
    s = PrintAmberTopFieldValues(lst,10,"%8i")
    string += s


#%FLAG NONBONDED_PARM_INDEX                                                      
#%FORMAT(10I8)                                                                   
#       1       2       4       7      11      16      22      29      37       2
#       3       5       8      12      17      23      30      38       4       5
#       6       9      13      18      24      31      39       7       8       9
#      10      14      19      25      32      40      11      12      13      14
#      15      20      26      33      41      16      17      18      19      20
#      21      27      34      42      22      23      24      25      26      27
#      28      35      43      29      30      31      32      33      34      35
#      36      44      37      38      39      40      41      42      43      44
#      45
    string+="%FLAG NONBONDED_PARM_INDEX\n"                                                      
    string+="%FORMAT(10I8)\n"                                                                   
    cnt=0
    M = [[0]*NTYPES for x in xrange(NTYPES)]
    for i in range(NTYPES):
        for j in range(i+1):
            cnt+=1
            iln = j
            icl = i
            M[iln][icl] = cnt
            M[icl][iln] = cnt
    M=np.array(M)
    lst=M.ravel() 
    s = PrintAmberTopFieldValues(lst,10,"%8i")
    string += s
        

#%FLAG RESIDUE_LABEL                                                             
#%FORMAT(20a4)                                                                   
#MOL 
    string+="%FLAG RESIDUE_LABEL\n"                                                             
    string+="%FORMAT(20a4)\n"                                                                   
    s = PrintAmberTopFieldValues([Mol.Name[0:4].upper()],20,"%4s")
    string += s

#%FLAG RESIDUE_POINTER                                                           
#%FORMAT(10I8)                                                                   
    string+="%FLAG RESIDUE_POINTER\n"
    string+="%FORMAT(10I8)\n"                                                                   
    s = PrintAmberTopFieldValues([1],10,"%8i")
    string += s

    
    if lBonded != "yes":
        s=""
        s+="%FLAG BOND_FORCE_CONSTANT\n" 
        s+="%FORMAT(5E16.8)\n\n"                                                    

        s+= "%FLAG BOND_EQUIL_VALUE\n"                                                       
        s+="%FORMAT(5E16.8)\n\n"                                                       

        s+="%FLAG ANGLE_FORCE_CONSTANT\n"                                                     
        s+="%FORMAT(5E16.8)\n\n"                                                      

        s+="%FLAG ANGLE_EQUIL_VALUE\n"                                                      
        s+="%FORMAT(5E16.8)\n\n"                                                      
        string+=s



    
    if lBonded != "yes":
        # To make it possible to calculate the 1-4 interactions we include the data about dihedral angles to topology. But
        # since no bonded interactions should be calculated, we ignore the functional type of the dihedrals in the Mol instance,
        # set the force constants to 0 for all dihedrals (get only the total number of torsions + impropers from the Mol instance)
        s=""
        s+="%FLAG DIHEDRAL_FORCE_CONSTANT\n"                                                   
        s+="%FORMAT(5E16.8)\n"                                                  
        lst = [0.0]*N_DIH_TYPES
        s += PrintAmberTopFieldValues(lst,5,"%16.8e")

        # Set periodicity to 1.0
        s+="%FLAG DIHEDRAL_PERIODICITY\n"                                                   
        s+="%FORMAT(5E16.8)\n"                                                   
        lst = [1.0]*N_DIH_TYPES
        s += PrintAmberTopFieldValues(lst,5,"%16.8e")

        # Set phases to 0.0
        s+="%FLAG DIHEDRAL_PHASE\n"                                                   
        s+="%FORMAT(5E16.8)\n"                                                   
        lst = [0.0]*N_DIH_TYPES
        s += PrintAmberTopFieldValues(lst,5,"%16.8e")
        
        string+=s



#%FLAG SCEE_SCALE_FACTOR                                                         
#%FORMAT(5E16.8)                                                                 
#  1.20000000E+00  1.20000000E+00  1.20000000E+00  1.20000000E+00  1.20000000E+00
#  1.20000000E+00  1.20000000E+00  1.20000000E+00  0.00000000E+00  0.00000000E+00
    string+="%FLAG SCEE_SCALE_FACTOR\n"                                                         
    string+="%FORMAT(5E16.8)\n"                                                                 
    lst = [1.0/FF.fudgeQQ]*N_DIH_TYPES
    s = PrintAmberTopFieldValues(lst,5,"%16.8e")
    string += s


#%FLAG SCNB_SCALE_FACTOR                                                         
#%FORMAT(5E16.8)                                                                 
#  2.00000000E+00  2.00000000E+00  2.00000000E+00  2.00000000E+00  2.00000000E+00
#  2.00000000E+00  2.00000000E+00  2.00000000E+00  0.00000000E+00  0.00000000E+00
    string+="%FLAG SCNB_SCALE_FACTOR\n"                                                         
    string+="%FORMAT(5E16.8)\n"                                                                 
    lst = [1.0/FF.fudgeLJ]*N_DIH_TYPES
    s = PrintAmberTopFieldValues(lst,5,"%16.8e")
    string += s


#%FLAG SOLTY                                                                     
#%FORMAT(5E16.8)                                                                 
#  0.00000000E+00  0.00000000E+00  0.00000000E+00  0.00000000E+00  0.00000000E+00
#  0.00000000E+00  0.00000000E+00  0.00000000E+00  0.00000000E+00  0.00000000E+00
    string+="%FLAG SOLTY\n"                                                                     
    string+="%FORMAT(5E16.8)\n"                                                                 
    lst = [0.0]*NTYPES
    s = PrintAmberTopFieldValues(lst,5,"%16.8e")
    string += s


#%FLAG LENNARD_JONES_ACOEF                                                       
#%FORMAT(5E16.8)                                                                 
#  3.79876399E+05  5.74393458E+05  8.19971662E+05  6.47841731E+05  9.24822270E+05
#  1.04308023E+06  5.44261042E+04  8.61541883E+04  9.71708117E+04  7.51607703E+03
#  6.06829342E+05  8.82619071E+05  9.95480466E+05  8.96776989E+04  9.44293233E+05
#  1.02595236E+03  2.27577561E+03  2.56678134E+03  1.07193646E+02  2.12601181E+03
#  1.39982777E-01  4.77908183E+04  7.62451550E+04  8.59947003E+04  6.55825601E+03
#  7.91627154E+04  8.90987508E+01  5.71629601E+03  4.71003287E+05  7.01803794E+05
#  7.91544157E+05  6.82786631E+04  7.44975864E+05  1.40467023E+03  6.00750218E+04
#  5.81803229E+05  0.00000000E+00  0.00000000E+00  0.00000000E+00  0.00000000E+00
#  0.00000000E+00  0.00000000E+00  0.00000000E+00  0.00000000E+00  0.00000000E+00
    string+="%FLAG LENNARD_JONES_ACOEF\n"                                                       
    if lChamber != 'yes':
        string+="%FORMAT(5E16.8)\n"                                                                 
    else:
        string+="%FORMAT(3E24.16)\n"                                                                 
        
    cnt=0
    Aij=[]
    Bij=[]
    Aij_14=[]
    Bij_14=[]
    for i in range(len(D_ATNUM_AT.keys())):
        for j in range(i+1):
            cnt+=1
            AT1 = FindObjByFieldValue( FF.AtomTypes, 'Type', D_ATNUM_AT[i+1] )
            AT2 = FindObjByFieldValue( FF.AtomTypes, 'Type', D_ATNUM_AT[j+1] )
            Func = CheckIfAllElementsInArrayAreTheSameAndReturn([AT1.Func,AT2.Func])
            p = ApplyCombRules( AT1.p, AT2.p, FF.CombRule, Func )
            p = ConvertAtomFuncParam(p,Func,'LJ_A_B')
            Aij.append(  p[0] )
            Bij.append(  p[1] )

            p_14 = ApplyCombRules( AT1.p_14, AT2.p_14, FF.CombRule, Func )
            p_14 = ConvertAtomFuncParam(p_14,Func,'LJ_A_B')
            Aij_14.append(  p_14[0] )
            Bij_14.append(  p_14[1] )

    lst = Aij  
    if lChamber != 'yes':
        s = PrintAmberTopFieldValues(lst,5,"%16.8e")
    else:
        s = PrintAmberTopFieldValues(lst,3,"%24.16e")
    string += s


#%FLAG LENNARD_JONES_BCOEF                                                       
#%FORMAT(5E16.8)                                                                 
#  5.64885984E+02  5.55666448E+02  5.31102864E+02  6.26720080E+02  5.99015525E+02
#  6.75612247E+02  1.11805549E+02  1.12529845E+02  1.26919150E+02  2.17257828E+01
#  6.77220874E+02  6.53361429E+02  7.36907417E+02  1.36131731E+02  8.01323529E+02
#  1.53505284E+01  1.82891803E+01  2.06278363E+01  2.59456373E+00  2.09604198E+01
#  9.37598976E-02  1.03580945E+02  1.04660679E+02  1.18043746E+02  2.00642027E+01
#  1.26451907E+02  2.33864085E+00  1.85196588E+01  6.29300710E+02  6.14502845E+02
#  6.93079947E+02  1.25287818E+02  7.50714425E+02  1.79702257E+01  1.16187983E+02
#  6.99746810E+02  0.00000000E+00  0.00000000E+00  0.00000000E+00  0.00000000E+00
#  0.00000000E+00  0.00000000E+00  0.00000000E+00  0.00000000E+00  0.00000000E+00
    string+="%FLAG LENNARD_JONES_BCOEF\n"                                                       
    if lChamber != 'yes':
        string+="%FORMAT(5E16.8)\n"                                                                 
    else:
        string+="%FORMAT(3E24.16)\n"                                                                 
    lst = Bij  
    if lChamber != 'yes':
        s = PrintAmberTopFieldValues(lst,5,"%16.8e")
    else:
        s = PrintAmberTopFieldValues(lst,3,"%24.16e")
    string += s

    if lChamber == 'yes':
        string+="%FLAG LENNARD_JONES_14_ACOEF\n"                                                       
        string+="%FORMAT(3E24.16)\n"
        lst = Aij_14  
        s = PrintAmberTopFieldValues(lst,3,"%24.16e")
        string += s

        string+="%FLAG LENNARD_JONES_14_BCOEF\n"                                                       
        string+="%FORMAT(3E24.16)\n"
        lst = Bij_14  
        s = PrintAmberTopFieldValues(lst,3,"%24.16e")
        string += s






    if lBonded != "yes":
        s=""
        s+="%FLAG BONDS_INC_HYDROGEN\n"                                                        
        s+="%FORMAT(10I8)\n\n"                                                                   

        s+="%FLAG BONDS_WITHOUT_HYDROGEN\n"                                                    
        s+="%FORMAT(10I8)\n\n"                                                                   

        s+="%FLAG ANGLES_INC_HYDROGEN\n"                                                       
        s+="%FORMAT(10I8)\n\n"                                                                   

        s+="%FLAG ANGLES_WITHOUT_HYDROGEN\n"                                                   
        s+="%FORMAT(10I8)\n\n"                                                                   


#%FLAG DIHEDRALS_INC_HYDROGEN                                                    
#%FORMAT(10I8)                                                                   
    s+="%FLAG DIHEDRALS_INC_HYDROGEN\n"                                                    
    s+="%FORMAT(10I8)\n"                                                                   
    lst=[]
    cnt=0
    used_dihedrals=[]
    used_pairs=[]

    for F in DIHEDRALS_H: 
        for FT in FF.TorsionTypes + FF.ImproperTypes:
                if (F.AtTypes == FT.AtTypes or F.AtTypes[::-1] == FT.AtTypes) and FT.__class__.__name__[:-4] == F.__class__.__name__:
                    cnt+=1
                    tmp = copy(F.AtNumbers)
                    AN = [ 3*(tmp[0]-1), 3*(tmp[1]-1), 3*(tmp[2]-1), 3*(tmp[3]-1) ]

                    if AN[2] == 0 or AN[3] == 0:
                        AN = AN[::-1]
                    
                    if tmp in used_dihedrals or tmp[::-1] in used_dihedrals:
                        AN[2] = - abs(AN[2])
                    else:
                        used_dihedrals.append(tmp)
                    
                    tmp = [ F.AtNumbers[0]  , F.AtNumbers[3]]
                    if tmp in used_pairs or tmp[::-1] in used_pairs:
                        AN[2] = - abs(AN[2])
                    else:
                        used_pairs.append(tmp)

                    if F.__class__.__name__ == 'Improper':
                        AN[3] = - abs(AN[3])
                        AN[2] = - abs(AN[2])

                    lst += [ AN[0], AN[1], AN[2], AN[3], cnt  ]
                    del(AN)
                    
    s += PrintAmberTopFieldValues(lst,10,"%8i")
    CNT_DA = cnt


#%FLAG DIHEDRALS_WITHOUT_HYDROGEN                                                
#%FORMAT(10I8)                                                                   
    s+="%FLAG DIHEDRALS_WITHOUT_HYDROGEN\n"                                                
    s+="%FORMAT(10I8)\n"                                                                   
    lst=[]
    cnt = CNT_DA
    #used_dihedrals=[]
    #used_pairs[]
    for F in DIHEDRALS_A: 
        for FT in FF.TorsionTypes + FF.ImproperTypes:
                if (F.AtTypes == FT.AtTypes or F.AtTypes[::-1] == FT.AtTypes) and FT.__class__.__name__[:-4] == F.__class__.__name__:
                #if F.AtTypes == FT.AtTypes or F.AtTypes[::-1] == FT.AtTypes:
                    cnt+=1
                    tmp = copy(F.AtNumbers)
                    AN = [ 3*(tmp[0]-1), 3*(tmp[1]-1), 3*(tmp[2]-1), 3*(tmp[3]-1) ]

                    if AN[2] == 0 or AN[3] == 0:
                        AN = AN[::-1]
                    
                    if tmp in used_dihedrals or tmp[::-1] in used_dihedrals:
                        AN[2] = - abs(AN[2])
                    else:
                        used_dihedrals.append(tmp)
                    
                    tmp = [ F.AtNumbers[0]  , F.AtNumbers[3]]
                    if tmp in used_pairs or tmp[::-1] in used_pairs:
                        AN[2] = - abs(AN[2])
                    else:
                        used_pairs.append(tmp)

                    if F.__class__.__name__ == 'Improper':
                        AN[3] = - abs(AN[3])
                        AN[2] = - abs(AN[2])

                    lst += [ AN[0], AN[1], AN[2], AN[3], cnt  ]
                    del(AN)
                    
    s += PrintAmberTopFieldValues(lst,10,"%8i")
    string += s



#%FLAG EXCLUDED_ATOMS_LIST
#%FORMAT(10I8)  (INB(i), i=1,NNB)
#  INB    : the excluded atom list.  To get the excluded list for atom 
#           "i" you need to traverse the NUMEX list, adding up all
#           the previous NUMEX values, since NUMEX(i) holds the number
#           of excluded atoms for atom "i", not the index into the 
#           NATEX list.  Let IEXCL = SUM(NUMEX(j), j=1,i-1), then
#           excluded atoms are INB(IEXCL) to INB(IEXCL+NUMEX(i)). Note,
#           this array was called NATEX at one point, and while in most
#           places it is now INB, it is still called NATEX in some places
#           (especially in pmemd)
#       2       3       4       5       6       7       8       9       3       4
#       5       6       7       8       9      10      19       4       5       6
#       7       8       9       5       6       7       6       7       7       8
#       9      10      11      12      17      19      20       9      10      19
#      10      11      12      13      14      17      18      19      20      11
#      12      13      14      15      17      19      20      12      13      14
#      19      13      14      15      16      17      18      19      14      15
#      17      15      16      17      18      19      20      16      17      18
#      19      17      18      19      20      19      20      20       0

#       2       3       4       5       6       7       8       9       3       4
#       5       6       7       8       9      10      19       4       5       6
#       7       8       9       5       6       7       6       7       7       8
#       9      10      11      12      17      19      20       9      10      19
#      10      11      12      13      14      17      18      19      20      11
#      12      13      14      15      17      19      20      12      13      14
#      19      13      14      15      16      17      18      19      14      15
#      17      15      16      17      18      19      20      16      17      18
#      19      17      18      19      20      19      20      20       0

    string+="%FLAG EXCLUDED_ATOMS_LIST\n"
    string+="%FORMAT(10I8)\n"
    
    lst=[]
    for k,v in Dexcl.iteritems():
        if len(v) !=0: 
            lst += sorted(v)
        else:
            lst += [ 0 ]
    s = PrintAmberTopFieldValues(lst,10,"%8i")
    string += s

#%FLAG HBOND_ACOEF                                                               
#%FORMAT(5E16.8)                                                                 
    string+="%FLAG HBOND_ACOEF\n"                                                        
    string+="%FORMAT(5E16.8)\n\n"                                                                   

#%FLAG HBOND_BCOEF                                                               
#%FORMAT(5E16.8)                                                                 
    string+="%FLAG HBOND_BCOEF\n"                                                               
    string+="%FORMAT(5E16.8)\n\n"                                                                 

#%FLAG HBCUT                                                                     
#%FORMAT(5E16.8)                                                                 
    string+="%FLAG HBCUT\n"                                                                     
    string+="%FORMAT(5E16.8)\n\n"                                                                 


#%FLAG AMBER_ATOM_TYPE                                                           
#%FORMAT(20a4)                                                                   
#o   c   c3  hc  hc  hc  n   hn  ca  ca  ha  ca  ha  ca  oh  ho  ca  ha  ca  ha  
    string+="%FLAG AMBER_ATOM_TYPE\n"                                                           
    string+="%FORMAT(20a4)\n"                                                                   
    lst=[]
    for A in Mol.Atoms: lst.append(A.Type[0:3].ljust(4))
    s = PrintAmberTopFieldValues(lst,20,"%4s")
    string += s

#%FLAG TREE_CHAIN_CLASSIFICATION                                                 
#%FORMAT(20a4)                                                                   
#M   M   3   E   E   E   M   E   M   B   E   B   E   B   S   E   S   E   M   E   

    string+="%FLAG TREE_CHAIN_CLASSIFICATION\n"                                                 
    string+="%FORMAT(20a4)\n"                                                                   
    lst = ['3']*NATOM
    s = PrintAmberTopFieldValues(lst,20,"%4s")
    string += s

#   try:
#       lst = Mol.ExtraData['TREE_CHAIN_CLASSIFICATION']
#       string+="%FLAG TREE_CHAIN_CLASSIFICATION\n"                                                 
#       string+="%FORMAT(20a4)\n"                                                                   
#       s = PrintAmberTopFieldValues(lst,20,"%4s")
#       string += s
#   except:
#       string+="%FLAG TREE_CHAIN_CLASSIFICATION\n"                                                 
#       string+="%FORMAT(20a4)\n\n"                                                                   
        


#%FLAG JOIN_ARRAY                                                                
#%FORMAT(10I8)                                                                   
#       0       0       0       0       0       0       0       0       0       0
#       0       0       0       0       0       0       0       0       0       0
    string+="%FLAG JOIN_ARRAY\n"                                                                
    string+="%FORMAT(10I8)\n"                                                                   
    lst = [0.0]*NATOM
    s = PrintAmberTopFieldValues(lst,10,"%8i")
    string += s



#%FLAG IROTAT                                                                    
#%FORMAT(10I8)                                                                   
#       0       0       0       0       0       0       0       0       0       0
#       0       0       0       0       0       0       0       0       0       0
    string+="%FLAG IROTAT\n"                                                                
    string+="%FORMAT(10I8)\n"                                                                   
    lst = [0.0]*NATOM
    s = PrintAmberTopFieldValues(lst,10,"%8i")
    string += s

#%FLAG RADIUS_SET                                                                
#%FORMAT(1a80)                                                                   
#modified Bondi radii (mbondi)                                                   
    string+="%FLAG RADIUS_SET\n"                                                                
    string+="%FORMAT(1a80)\n"                                                                   
    lst=['not given']
    s = PrintAmberTopFieldValues(lst,1,"%80s")
    string += s

#   try:
#       lst = Mol.ExtraData['BORN_RADIUS_SET']
#       string+="%FLAG RADIUS_SET\n"                                                                
#       string+="%FORMAT(1a80)\n"                                                                   
#       s = PrintAmberTopFieldValues(lst,1,"%80s")
#       string += s
#   except:
#       string+="%FLAG RADIUS_SET\n"                                                                
#       string+="%FORMAT(1a80)\n\n"                                                                   
        

#%FLAG RADII                                                                     
#%FORMAT(5E16.8)                                                                 
#  1.50000000E+00  1.70000000E+00  1.70000000E+00  1.30000000E+00  1.30000000E+00
#  1.30000000E+00  1.55000000E+00  1.30000000E+00  1.70000000E+00  1.70000000E+00
#  1.30000000E+00  1.70000000E+00  1.30000000E+00  1.70000000E+00  1.50000000E+00
#  8.00000000E-01  1.70000000E+00  1.30000000E+00  1.70000000E+00  1.30000000E+00
    string+="%FLAG RADII\n"                                                                
    string+="%FORMAT(5E16.8)\n"                                                                   
    lst=[0.0]*NATOM
    s = PrintAmberTopFieldValues(lst,5,"%16.8e")
    string += s
#   try:
#       lst = Mol.ExtraData['BORN_RADII']
#       string+="%FLAG RADII\n"                                                                
#       string+="%FORMAT(5E16.8)\n"                                                                   
#       s = PrintAmberTopFieldValues(lst,5,"%16.8e")
#       string += s
#   except:
#       string+="%FLAG RADII\n"                                                                
#       string+="%FORMAT(5E16.8)\n\n"                                                                   
        


#%FLAG SCREEN                                                                    
#%FORMAT(5E16.8)                                                                 
#  8.50000000E-01  7.20000000E-01  7.20000000E-01  8.50000000E-01  8.50000000E-01
#  8.50000000E-01  7.90000000E-01  8.50000000E-01  7.20000000E-01  7.20000000E-01
#  8.50000000E-01  7.20000000E-01  8.50000000E-01  7.20000000E-01  8.50000000E-01
#  8.50000000E-01  7.20000000E-01  8.50000000E-01  7.20000000E-01  8.50000000E-01
    string+="%FLAG SCREEN\n"                                                                
    string+="%FORMAT(5E16.8)\n"                                                                   
    lst=[0.0]*NATOM
    s = PrintAmberTopFieldValues(lst,5,"%16.8e")
    string += s

#   try:
#       lst = Mol.ExtraData['BORN_SCREEN']
#       string+="%FLAG SCREEN\n"                                                                
#       string+="%FORMAT(5E16.8)\n"                                                                   
#       s = PrintAmberTopFieldValues(lst,5,"%16.8e")
#       string += s
#   except:
#       string+="%FLAG SCREEN\n"                                                                
#       string+="%FORMAT(5E16.8)\n\n"                                                                   
        

    of.write(string)















