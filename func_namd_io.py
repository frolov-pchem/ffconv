from molecule_class import *
from copy import copy

NaN=float('NaN')


### TODO: Old function: not tested recently, might not work due to some changes in molecule_class.py
def ReadNamdPSFFile(fn):
    
    AtomList = [] ; AtType=[] ;  AtName=[] ; q = [] ; Mass = [] ; ResidueName=[]; ResidueNumber=[]; MolName=[]
    Pair14List = []
    BondList = [] ; BdNumbers = [] 
    AngleList = [] ; AnNumbers = [] 
    TorsionList = [] ; TnNumbers = [] 
    ImproperList = [] ; ImNumbers = [] 

    FieldsList=['!NATOM','!NBOND:','!NTHETA:','!NPHI:','!NIMPHI:']

    ifile=open(fn,'r')
    iBond = 0 ; iAngle = 0 ; iTorsion = 0 ; iImproper=0
    #lines=ifile.readlines()
    while 1:
        line=ifile.readline()
    
    #for line in lines:
        if not line: break
        ln=line.split()
        if len(ln)<2: continue
        
        if ln[1] in FieldsList:
            field=ln[1]
        else:
            field='none'
        
#      33 !NATOM
#       1 IBPR 1    IBPR O1   O11   -0.636249       15.9940           0   0.00000     -0.301140E-02
#       2 IBPR 1    IBPR C1   CA1    0.226021       12.0110           0   0.00000     -0.301140E-02
        if field=='!NATOM':
            NAtoms=int(ln[0])
            print NAtoms
            for i in range(NAtoms):
                line=ifile.readline()
                ln=line.split()
                AtomList.append(int(ln[0]))
                ResidueNumber.append(int(ln[2]))
                ResidueName.append(ln[3])
                MolName.append(ln[1])
                AtName.append(ln[4])
                AtType.append(ln[5])
                q.append(float(ln[6]))
                Mass.append(float(ln[7]))
       
#      33 !NBOND: bonds
#       1      10       1      16       2       7       2       8
#       3       4       3      17       4      18       5       6
        elif field=='!NBOND:':
            NBonds=int(ln[0])
            line=ifile.readline(); ln=line.split()
            lns=ln
            while len(ln)!=0:
                line=ifile.readline(); ln=line.split()
                lns+=ln 
            for i in range(NBonds):
                iBond += 1
                BondList.append(iBond)
                BdNumbers.append( [ int(lns[2*i]),int(lns[2*i+1]) ] )

#     116 !NTHETA: angles
#      10       1      16      10       1      16       3       2       7
#       3       2       7       3       2       8       3       2       8
        elif field=='!NTHETA:':
            NAngles=int(ln[0])
            line=ifile.readline(); ln=line.split()
            lns=ln
            while len(ln)!=0:
                line=ifile.readline(); ln=line.split()
                lns+=ln 
            for i in range(NAngles):
                iAngle += 1
                AngleList.append(iAngle)
                AnNumbers.append( [ int(lns[3*i]),int(lns[3*i+1]),int(lns[3*i+2]) ] )

#      80 !NPHI: dihedrals
#       1      10       9       5       1      10       9      13
#       1      10       9      23       2       3       4       5
        elif field=='!NPHI:':
            NTorsions=int(ln[0])
            line=ifile.readline(); ln=line.split()
            lns=ln
            while len(ln)!=0:
                line=ifile.readline(); ln=line.split()
                lns+=ln 
            for i in range(NTorsions):
                iTorsion += 1
                TorsionList.append(iTorsion)
                TnNumbers.append( [ int(lns[4*i]),int(lns[4*i+1]),int(lns[4*i+2]),int(lns[4*i+3]) ] )

#        1 !NIMPHI: impropers
#       10         9        11         1

        elif field=='!NIMPHI:':
            NImpropers=int(ln[0])
            line=ifile.readline(); ln=line.split()
            lns=ln
            while len(ln)!=0:
                line=ifile.readline(); ln=line.split()
                lns+=ln 
            for i in range(NImpropers):
                iImproper += 1
                ImproperList.append(iImproper)
                ImNumbers.append( [ int(lns[4*i]),int(lns[4*i+1]),int(lns[4*i+2]),int(lns[4*i+3]) ] )

        else:
            continue

    ifile.close()
    
    #print ResidueName, ResidueNumber, MolName
    #print TnNumbers
    #print ImNumbers
    #ImNumbers = []
    ResName = CheckIfAllElementsInArrayAreTheSameAndReturn(ResidueName)
    ResNumber = CheckIfAllElementsInArrayAreTheSameAndReturn(ResidueNumber)
    MoleculeName = CheckIfAllElementsInArrayAreTheSameAndReturn(MolName)
    
    #print ResName, ResNumber, MoleculeName

    if ResName != 'none' and ResNumber != 'none' and MoleculeName != 'none':
        Mol = Molecule(MoleculeName)
    else:
        print "!!! Error: There are several residue names, residue numbers or molecule names in the PSF file. The script supports only single molecules. Exiting."
        return 'none'
    Pr14Numbers=[]
    Mol.FillMoleculeTopology(AtType,AtName,q,BdNumbers,Pr14Numbers,AnNumbers,TnNumbers,ImNumbers)
    
    #print Mol
    return Mol









