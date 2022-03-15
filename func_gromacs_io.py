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
import sys
import datetime
import getpass

# TODO: extend the section fileds headings for each section name [ section ]
# - finish ReadGromacsMoleculeITPFile    
# - print different heading depending on FuncType. Add units to heading.
# - Check for NaN output ?

###
### Set global variables
### 
user=str(getpass.getuser())
#today=str(datetime.date.today())
#today = datetime.today()
today=str(today.strftime('%d %b %Y, %H:%M:%S'))


### Gromacs specific global variables

# Known Gromacs interaction function types (consult Gromacs Manual).  
# TODO: extend the list according to GMX manual!
GMXNBFuncTypes = {'Default':1, 'LJ_sig_eps':1, 'Buckingham':2}
GMXCombRuleTypes = {'Default':3, 'LorentzBerthelot':2 , 'Geometric':3, 'GeometricAB':1 }
GMXBondFuncTypes = {'Default':1, 'Harmonic':1} # To be extended
GMXAngleFuncTypes = {'Default':1, 'Harmonic':1, 'UreyBradley':5} # To be extended
GMXTorsionFuncTypes = {'Default':1, 'Periodic':1, 'PeriodicMultiple':9, 'Fourier':5, 'RB':3} # To be extended
GMXImproperFuncTypes = {'Default':2, 'Periodic':4, 'PeriodicMultiple':9, 'Harmonic':2, 'Fourier':5, 'RB':3} # To be extended

GMXFuncTypes={}
GMXFuncTypes['NB'] = GMXNBFuncTypes
GMXFuncTypes['CombRule'] = GMXCombRuleTypes
GMXFuncTypes['Bond'] = GMXBondFuncTypes
GMXFuncTypes['Angle'] = GMXAngleFuncTypes
GMXFuncTypes['Torsion'] = GMXTorsionFuncTypes
GMXFuncTypes['Improper'] = GMXImproperFuncTypes





#----------------------------------------------------------------------------------------------------------------------
def WriteGromacsFFFile(FName,FF,**kwargs):
#
#   EXAMPLE USAGE
# WriteGromacsFFFile('paracetamol_AT.itp',FF,lBonded='yes',NBFuncType='LJ_sig_eps',lGenPairs='yes',lDef='no')
#
#   DESCRIPTION
# The function writes the include force field file in Gromacs format (itp file) from a ForceField object. 
# The following fields will be created if present in the data present in FF object:
# [ defaults ], [ atomtypes ], [ pairtypes ], [ bondtypes], [ dihedrals ], [ dihedrals ] ; improper
# 
# if lBonded != 'yes':
#       only [ defaults ] and [ atomtypes ] sections will be printed. This can be convinient when in a molecular itp
#       file the potential parameters are given explicitly. (Mind! only Atom parameters could not be set explicitly 
#       in the molecule itp file).
#
# if lDef == 'yes':
#       the [ defaults ] section will be commented out.
#
#   
#   INPUT
# Fname - name of the output file (string)
# FF - a ForceField class object
# (Optional)
# lBonded - (default 'yes') Key to write the bonded potentials ('yes'), or skip them (print only atomtypes) ('no')
# lDef - (default 'no') Key to write the [ defaults ] force field section ('yes'), or write it commentrd out ('no')
# lGenPairs - (default depends on FF.Family) Key to write the third item of the [ defaults ] force field section. ('yes') - ask to generate pair potential parameters, or not to - ('no').
# NBFuncType - (default depends on FF.Family) the functional form for nonbonded interactions (used to convert the FF data to the required form) (string)
#              Must be one of the values in the global dictinary GMXNBFuncTypes. 
#
#   REQUIREMENTS
# The FF object AtomTypes[i].Func must be the same (otherwise error returned).
#
    ### TODO: 
        # Not fully tested
        # Write documentation
    
    ###
    ### Necessary input checks
    ###
    if FF.__class__.__name__ != 'ForceField':
            sys.stderr.write('--- !Error in function '+inspect.stack()[0][3]+': the Second function argument must be of class ForceField. Class of the given object: ['+FF.__class__.__name__+']. Aborting.\n')
            sys.exit(1)

    if FName.__class__.__name__ != 'str':
            sys.stderr.write('--- !Error in function '+inspect.stack()[0][3]+': the First function argument must be of class str. Class of the given object: ['+FName.__class__.__name__+']. Aborting.\n')
            sys.exit(1)


    ###
    ### Defaults
    ### 
    Vars={}
    Vars['lBonded']='yes'
    Vars['lDef']='no'
    Vars['ofdef']='defaults.itp'

    
    # Set default values for vars:
    if   FF.Family.upper() == 'CHARMM':
        Vars['NBFuncType']='LJ_sig_eps' ; Vars['lGenPairs']='yes'
    elif FF.Family.upper() == 'OPLS':
        Vars['NBFuncType']='LJ_sig_eps' ; Vars['lGenPairs']='yes'
    elif FF.Family.upper() == 'AMBER':
        Vars['NBFuncType']='LJ_sig_eps' ; Vars['lGenPairs']='yes'
    elif FF.Family.upper() == 'none':
        Vars['NBFuncType']='LJ_sig_eps' ; Vars['lGenPairs']='yes'
        sys.stderr.write('--- !Warning in function '+inspect.stack()[0][3]+': FF.Family ['+FF.Family+'] is set to none. Set NBFuncType=LJ_sig_eps. lGenPairs=yes Continue.')
    else:
        Vars['NBFuncType']='LJ_sig_eps' ; Vars['lGenPairs']='yes'
        sys.stderr.write('--- !Warning in function '+inspect.stack()[0][3]+': FF.Family ['+FF.Family+'] is NOT known. Set NBFuncType=LJ_sig_eps, lGenPairs=yes. Continue.')

    ###
    ### Process function arguments 
    ###
    for k,v in kwargs.iteritems():
        if not k in Vars.keys():
            sys.stderr.write('--- !Warning in function '+inspect.stack()[0][3]+' not known input key argument '+k+'. Continue.\n')
        for key,val in Vars.iteritems():
            if k==key:
                Vars[k]=v
            else:
                pass

    ###
    ### Set the variables
    ###
    NBFuncType = Vars['NBFuncType']
    lBonded = Vars['lBonded']
    lGenPairs = Vars['lGenPairs']
    lDef = Vars['lDef']
    ofdef = Vars['ofdef']

    # Set FF variables
    fudgeLJ = FF.fudgeLJ
    fudgeQQ = FF.fudgeQQ
    #NREXCL = FF.NREXCL  # not used in this function
    CombRule = FF.CombRule
    

    ###
    ### Checks
    ###
    if lGenPairs != 'yes' and lGenPairs != 'no':
            sys.stderr.write('--- !Error in function '+inspect.stack()[0][3]+': the lGenPairs must be yes or no. ['+lGenPairs+'] given. Aborting.\n')
            sys.exit(1)
        
    if lBonded != 'yes' and lBonded != 'no':
            sys.stderr.write('--- !Error in function '+inspect.stack()[0][3]+': the lBonded must be yes or no. ['+lBonded+'] given. Aborting.\n')
            sys.exit(1)
        
    if lDef != 'yes' and lDef != 'no':
            sys.stderr.write('--- !Error in function '+inspect.stack()[0][3]+': the lDef must be yes or no. ['+lDef+'] given. Aborting.\n')
            sys.exit(1)
        
    if not CombRule in GMXCombRuleTypes.keys():
            sys.stderr.write('--- !Warning in function '+inspect.stack()[0][3]+': CombRule ['+CombRule+'] is not present in GMXCombRuleTypes global dictionary. CombRule set to NaN. Continue\n')
            CombRule=float('NaN')
        
    if not NBFuncType in GMXNBFuncTypes.keys():
            sys.stderr.write('--- !Error in function '+inspect.stack()[0][3]+': the NBFuncType must be present in the GMXNBFuncTypes global dictionary. ['+NBFuncType+'] given, the dictionary GMXNBFuncTypes ['+str(GMXNBFuncTypes)+']. Extend the dictionary if necessary. Aborting.\n')
            sys.exit(1)
    
    if fudgeLJ < 0 or fudgeLJ > 1:
            sys.stderr.write('--- !Warning in function '+inspect.stack()[0][3]+': fudgeLJ is usually set in the interval [0:1], but ['+str(fudgeLJ)+'] given. Check the data. Continue\n')
    if fudgeQQ < 0 or fudgeQQ > 1:
            sys.stderr.write('--- !Warning in function '+inspect.stack()[0][3]+': fudgeQQ is usually set in the interval [0:1], but ['+str(fudgeQQ)+'] given. Check the data. Continue\n')
        

    ###
    ### Perform conversions
    ###
    
    # Convert. This also takes care that all nonbonded functional types are the same
    ConvertNonBondedParamsTo(FF,NBFuncType)
    
    # Convert units
    ConvertUnitsTo(FF,'togmx')
    
    # Convert angle units to degree 
    ConvertUnitsTo(FF,'togmx',lrad2degree='yes')


    # This required to trasform e.g. CHARMM shortcats to explicit names
    FF.ChangeJointAtomTypeToAtomTypesBelongingToIt() # change atomtypes with '*' to the atom types belonging to it, eg 'CA*' to 'CA1, CA2, ...'


    ###
    ### Writing the FF itp
    ###
    
##define _FF_CHARMM
#[ defaults ]
#; nbfunc	comb-rule	gen-pairs	fudgeLJ	fudgeQQ
#1	2	yes	1.0	1.0
    of=open(FName,'w')

    # Header
    s=";;; This is the itp for ForceField ["+FF.Name+"] created by the ["+sys.argv[0]+"] tool "
    s2=";;; Generated by user " +user+ " at "+today
    of.write(";"*len(s)+"\n"\
            +s+"\n"\
            +s2+"\n"\
            +";"*len(s)+"\n\n")


    # Writing FFCONV info
    string= \
    "; -- ffconv.py info: \n" 
    for k,v in Vars.iteritems():
        string += "; "+str(k)+" = "+str(v)+"\n"
    string+= "\n"
    
    # Writing FF info
    string+= \
    "; -- force field parameters: \n" +\
    "; Name = "+str(FF.Name)+"\n" +\
    "; Family = "+str(FF.Family)+"\n" +\
    "; Units = "+str(FF.Units)+"\n" +\
    "; GeneratedBy = "+str(FF.GeneratedBy)+"\n" +\
    "; Comment = "+str(FF.Comment)+"\n" +\
    "; NREXCL = "+str(FF.NREXCL)+"\n" +\
    "; fudgeLJ = "+str(FF.fudgeLJ)+"\n" +\
    "; fudgeQQ = "+str(FF.fudgeQQ)+"\n" +\
    "; CombRule = "+str(FF.CombRule)+"\n" +\
    "\n" 

    # Force field defaults section
#   if lDef == 'no':
#       string+= \
#   ";;; compatible [ defaults ] section for gromacs topology. \n" +\
#   ";[ defaults ]\n" +\
#   ";; nbfunc comb-rule  gen-pairs   fudgeLJ   fudgeQQ\n" +\
#   "; "+str(GMXNBFuncTypes[NBFuncType])+"    "+str(GMXCombRuleTypes[CombRule])+"    "+str(lGenPairs)+"    "+str(fudgeLJ)+"   "+str(fudgeQQ)+"\n\n\n" 
#   else:
#       string+= \
#   "[ defaults ]\n" +\
#   "; nbfunc comb-rule  gen-pairs   fudgeLJ   fudgeQQ\n" +\
#   " "+str(GMXNBFuncTypes[NBFuncType])+"    "+str(GMXCombRuleTypes[CombRule])+"    "+str(lGenPairs)+"    "+str(fudgeLJ)+"   "+str(fudgeQQ)+"\n\n\n" 

    string+= \
   ";;; compatible [ defaults ] section for this gromacs FF itp file. \n" +\
   ";[ defaults ]\n" +\
   ";; nbfunc comb-rule  gen-pairs   fudgeLJ   fudgeQQ\n" +\
   "; "+str(GMXNBFuncTypes[NBFuncType])+"    "+str(GMXCombRuleTypes[CombRule])+"    "+str(lGenPairs)+"    "+str(fudgeLJ)+"   "+str(fudgeQQ)+"\n\n\n" 

    
    if lDef == 'yes':
        # Print the FF defaults itp file
        defstring = \
   "[ defaults ]\n" +\
   "; nbfunc comb-rule  gen-pairs   fudgeLJ   fudgeQQ\n" +\
   " "+str(GMXNBFuncTypes[NBFuncType])+"    "+str(GMXCombRuleTypes[CombRule])+"    "+str(lGenPairs)+"    "+str(fudgeLJ)+"   "+str(fudgeQQ)+"\n\n\n" 
        ofdefaults=open(ofdef,'w')
        ofdefaults.write(defstring)

        

    


    # FF atom types
#[ atomtypes ]
#;name	at.num	mass	charge	ptype	sigma	epsilon
#C	6	12.01100	0.51	A	0.356359487256	0.46024 
#CA	6	12.01100	-0.115	A	0.355005321205	0.29288 
#CC	6	12.01100	0.62	A	0.356359487256	0.29288 
#CD	6	12.01100	0.000	A	0.356359487256	0.29288 ; partial charge def not found
    
    Fs=FF.AtomTypes
    if len(Fs) > 0:
        string += \
        "[ atomtypes ]\n"+\
        ";name	at.num	mass	charge	ptype	sigma[nm]	epsilon[kJ_mol]\n"
        for F in Fs:
            s="%6s %6i %10.4f %14.8f %4s" % (F.Type,F.PeriodicTableNum,F.Mass,0.0,'A') 
            for p in F.p:
                tmp="%14.8f" % (p)
                s += tmp
            if len(F.Description) != 0: s += "  ; desc: "+F.Description+""
            s+='\n'
            string+=s
        string+='\n'
    else:   
        sys.stderr.write('--- !Warning in function '+inspect.stack()[0][3]+': there are no AtomTypes for the given FF named ['+FF.Name+']. Continue\n')
        

    if lBonded == 'yes':
    
#[ pairtypes ]
#; i j   func    sigma1-4    epsilon1-4 ; THESE ARE 1-4 INTERACTIONS
#CP1 CP1 1   0.338541512893  0.04184
#CP1 CP2 1   0.338541512893  0.04184
      Fs=FF.PairTypes
      if len(Fs) > 0:
        string += \
        "[ pairtypes ]\n"+\
        "; i    j    func     sig1-4    eps1-4 [the line is more likely correct]\n"
        for F in Fs:
                 
            s="%6s %6s %4i" % (F.AtTypes[0],F.AtTypes[1],GMXNBFuncTypes[NBFuncType]) 
            for p in F.p:
                tmp="%16.6f" % (p)
                s += tmp
            if len(F.Description) != 0: s += "  ; desc: "+F.Description+""
            s+='\n'
            string+=s
        string+='\n'





#[ bondtypes ]
#; i	j	func	b0	kb
#CST	OST	1	0.116	784884.928
#SS	FE	1	0.232	209200.0

      Fs=FF.BondTypes
      if len(Fs) > 0:
        string += \
        "[ bondtypes ]\n"+\
        "; i    j    func     b0    kb [the line is more likely correct]\n"
        for F in Fs:
            if not F.Func in GMXFuncTypes['Bond'].keys():
                sys.stderr.write('--- !Error in function '+inspect.stack()[0][3]+': the Bond funtion type ['+F.Func+'] is not present in the GMXFuncTypes[Bond] global dictionary.\n The dictionary GMXFuncTypes[Bond] ['+str(GMXFuncTypes['Bond'])+']\n. BondType: ['+str(F.__dict__)+']\n. Extend the dictionary if necessary. Aborting.\n')
                sys.exit(1)
            s="%6s %6s %4i" % (F.AtTypes[0],F.AtTypes[1],GMXBondFuncTypes[F.Func]) 
            for p in F.p:
                tmp="%20.6f" % (p)
                s += tmp
            if len(F.Description) != 0: s += "  ; desc: "+F.Description+""
            s+='\n'
            string+=s
        string+='\n'

#[ angletypes ]
#; i	j	k	func	th0	cth	ub0	cub
#OST	CST	OST	5	180.0000	25104.0	0.0	0.0
#CS	SS	FE	5	100.6	418.4	0.0	0.0
      Fs=FF.AngleTypes
      if len(Fs) > 0:
        string += \
        "[ angletypes ]\n"+\
        "; i	j	k	func	th0	cth	ub0	cub [the line is more likely correct]\n" # In principle this should change depending on AN.Func. probably to be extended

        for F in Fs:
            if not F.Func in GMXFuncTypes['Angle'].keys():
                sys.stderr.write('--- !Error in function '+inspect.stack()[0][3]+': the Angle funtion type ['+F.Func+'] is not present in the GMXFuncTypes[Angle] global dictionary.\n The dictionary GMXFuncTypes[Angle] ['+str(GMXFuncTypes['Angle'])+']\n. AngleType: ['+str(F.__dict__)+']\n. Extend the dictionary if necessary. Aborting.\n')
                sys.exit(1)
            s="%6s %6s %6s %4i" % (F.AtTypes[0],F.AtTypes[1],F.AtTypes[2],GMXAngleFuncTypes[F.Func]) 
            for p in F.p:
                tmp="%16.6f" % (p)
                s += tmp
            if len(F.Description) != 0: s += "  ; desc: "+F.Description+""
            s+='\n'
            string+=s
        string+='\n'


#[ dihedraltypes ] ; proper dihedrals
#; i j   k   l   func    phi0    cp  mult
#C   CT1 NH1 C   9   180.00  0.8368  1
#C   CT2 NH1 C   9   180.00  0.8368  1
      Fs=FF.TorsionTypes
      if len(Fs) > 0:
        string += \
        "[ dihedraltypes ] ; proper dihedrals \n"+\
        "; i j   k   l   func    phi0    cp  mult [the line is more likely correct]\n" # In principle this should change depending on AN.Func. probably to be extended

        for F in Fs:
            if not F.Func in GMXFuncTypes['Torsion'].keys():
                sys.stderr.write('--- !Error in function '+inspect.stack()[0][3]+': the Torsion funtion type ['+F.Func+'] is not present in the GMXFuncTypes[Torsion] global dictionary.\n The dictionary GMXFuncTypes[Torsion] ['+str(GMXFuncTypes['Torsion'])+']\n. TorsionType: ['+str(F.__dict__)+']\n. Extend the dictionary if necessary. Aborting.\n')
                sys.exit(1)
            s="%6s %6s %6s %6s %4i" % (F.AtTypes[0],F.AtTypes[1],F.AtTypes[2],F.AtTypes[3],GMXTorsionFuncTypes[F.Func]) 
            for p in F.p:
                tmp="%16.6f" % (p)
                s += tmp
            if len(F.Description) != 0: s += "  ; desc: "+F.Description+""
            s+='\n'
            string+=s
        string+='\n'

#[ dihedraltypes ]
#; i j   k   l   func    q0  cq
#CPB CPA NPH CPA 2   0.0000  174.0544
#HA  CPA CPA CPM 2   0.0000  246.0192
      Fs=FF.ImproperTypes
      if len(Fs) > 0:
        string += \
        "[ dihedraltypes ] ; improper dihedrals\n"+\
        "; i   j   k   l   func    q0  cq [the line is more likely correct]\n" # In principle this should change depending on IM.Func. probably to be extended

        for F in Fs:
            if not F.Func in GMXFuncTypes['Improper'].keys():
                sys.stderr.write('--- !Error in function '+inspect.stack()[0][3]+': the Improper funtion type ['+F.Func+'] is not present in the GMXFuncTypes[Improper] global dictionary.\n The dictionary GMXFuncTypes[Improper] ['+str(GMXFuncTypes['Improper'])+']\n. ImproperType: ['+str(F.__dict__)+']\n. Extend the dictionary if necessary. Aborting.\n')
                sys.exit(1)
            s="%6s %6s %6s %6s %4i" % (F.AtTypes[0],F.AtTypes[1],F.AtTypes[2],F.AtTypes[3],GMXImproperFuncTypes[F.Func]) 
            for p in F.p:
                tmp="%16.6f" % (p)
                s += tmp
            if len(F.Description) != 0: s += "  ; desc: "+F.Description+""
            s+='\n'
            string+=s
        string+='\n'

    of.write(string)







#----------------------------------------------------------------------------------------------------------------------
def WriteGromacsMoleculeITPFile(FName,Mol,**kwargs):
#
#   EXAMPLE USAGE
# WriteGromacsMoleculeITPFile('paracetamol.itp',Mol,Style='Explicit',NREXCL=3,prfunc=1,bdfunc=1,anfunc=1,tnfunc=1,imfunc=1,NBFuncType='LJ_sig_eps')
#
#   DESCRIPTION
# The function writes the molecular include topology file in Gromacs format (itp file) from a Molecule object. 
#
# if Style = Implicit:
#    no potential parameter will be printed (in this case they are supposed to be taken by Gromacs 
#    from the Force Field ITP file based on the atom types of molecule atoms occuring in a bobded 
#    interaction). However, unfortunately, Gromacs requires the functional form of the all interactions to be printed
#    explicitly, anyway. 
#
#   if functional type(s) of the bonded interactions not given as arguments to function (e.g. bdfunc,anfunc,
#   tnfunc,imfunc = 'none'): 
#       The function will get the functional form from each entry of the bonded interactions from the 
#       Mol data (this must be in the GMXFuncTypes global dictionary). If a functional form for a bonded interaction is 
#       absent in Mol, then the script will put the default value given for each type of interaction in the global 
#       dictionary GMXFuncTypes (If a functional form for a pair interaction is absent in Mol, then the script will 
#       put the atom functional type given in Mol).
#   else:
#      The functional type(s) will be set the given value and the data of Mol ignored.
#   
# if Style = Explicit:
#    potential parameter will be printed to output.
#
#   if functional type(s) of the bonded interactions not given as arguments to function (e.g. bdfunc,anfunc,
#   tnfunc,imfunc = 'none'): 
#       The function will get the functional form from each entry of the bonded interactions from the 
#       Mol data (this must be in the GMXFuncTypes global dictionary). If a functional form for a bonded interaction is 
#       absent in Mol, then the script will terminate with error. (If a functional form for a pair interaction is absent
#       in Mol, then the script will put the atom functional type given in Mol if pair potential parametes are NaN, them
#       they are ignored in the output).
#   else:
#      The functional type(s) will be set the given value and the data of Mol ignored. The potential will be printed 
#      as they were.
#   
#   INPUT
# Fname - name of the output file (string)
# Mol - a Molecule class object
# (Optional)
# prfunc - (default 'none') functional form for pair interactions in Gromacs notations (integer)
# bdfunc - (default 'none') functional form for bond interactions in Gromacs notations (integer)
# anfunc - (default 'none') functional form for angle interactions in Gromacs notations (integer)
# tnfunc - (default 'none') functional form for torsion angle interactions in Gromacs notations (integer)
# imfunc - (default 'none') functional form for improper dihedral angle interactions in Gromacs notations (integer)
# NREXCL - (default depends on Mol.FFName) number to exclude nonbonded interaction between closely bonded atoms in molecule (integer)
#          e.g. 3 would mean that the interaction between atoms which are separated within one molecule no more than by 3 covalent bonds are excluded  
# Style - (default 'Implicit') Key to write the bonded interactions parameters potential parametes explicitly ('Explicit'), or omit them ('Implicit')
# NBFuncType - (default depends on Mol.FFName) the functional form for nonbonded interactions (used to convert the Mol data to the required form) (string)
#              Must be one of the values in the global dictinary GMXNBFuncTypes. 
#
#   REQUIREMENTS
# The Molecule object is supposed to have Atoms with Type, Name, Number, Charge, Mass specified attributes. 
# The Atoms function types must be the same (otherwise error returned).
#
#
    ### TODO: 
        # Not fully tested
    
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
    #Vars['NREXCL']=float('NaN')
    #Vars['NBFuncType']='none'
    Vars['Style']='Implicit'
    
    Vars['atfunc']='none'
    Vars['bdfunc']='none'
    Vars['anfunc']='none'
    Vars['tnfunc']='none'
    Vars['imfunc']='none'
    Vars['prfunc']='none'    


    # TODO rewrite with FFFamily
    # Set default values for Vars:

    FFFamily = GetFFFamily(Mol.FFName)
    if   FFFamily == 'CHARMM':
        Vars['NREXCL']=3 ;  Vars['NBFuncType']='LJ_sig_eps' 
    elif FFFamily == 'OPLS':
        Vars['NREXCL']=3 ;  Vars['NBFuncType']='LJ_sig_eps' 
    elif FFFamily == 'AMBER':
        Vars['NREXCL']=3 ;  Vars['NBFuncType']='LJ_sig_eps' 
    elif FFFamily == 'none':
        sys.stderr.write('--- !Warning in function '+inspect.stack()[0][3]+': FF family is none for Mol.FFName ['+Mol.FFName+']. Set NREXCL=3, NBFuncType=LJ_sig_eps. Continue.\n')
        Vars['NREXCL']=3 ;  Vars['NBFuncType']='LJ_sig_eps' 
    else:
        sys.stderr.write('--- !Warning in function '+inspect.stack()[0][3]+': FF family ['+FFFamily+'] is NOT known by the function. Set NREXCL=3, NBFuncType=LJ_sig_eps. Continue.\n')
        Vars['NREXCL']=3 ;  Vars['NBFuncType']='LJ_sig_eps' 


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

    NBFuncType=Vars['NBFuncType']
    NREXCL=int(float(Vars['NREXCL']))
    Style=Vars['Style']


    ###
    ### Checks
    ###
    if not NBFuncType in GMXNBFuncTypes.keys():
            sys.stderr.write('--- !Error in function '+inspect.stack()[0][3]+': the NBFuncType must be present in the GMXNBFuncTypes global dictionary. ['+NBFuncType+'] given, the dictionary GMXNBFuncTypes ['+str(GMXNBFuncTypes)+']. Extend the dictionary if necessary. Aborting.\n')
            sys.exit(1)
        
        
    if Style != 'Explicit' and Style != 'Implicit':
            sys.stderr.write('--- !Error in function '+inspect.stack()[0][3]+': the Style to write molecular ITP file must be Implicit or Explicit, but ['+Style+'] given. The molecule instance named ['+Mol.Name+']. Aborting.\n')
            sys.exit(1)
        
    if NREXCL < 0 or NREXCL > 3:
            sys.stderr.write('--- !Warning in function '+inspect.stack()[0][3]+': NREXCL is usually set in the interval [0:3], but ['+str(NREXCL)+'] given. Check the data. Continue\n')
        

    if len(Mol.Atoms) < 1:
            sys.stderr.write('--- !Error in function '+inspect.stack()[0][3]+': the number of atoms for the molecule named ['+Mol.Name+'] = 0. Not allowed. Aborting.\n')
            sys.exit(1)
        


        
    
    ###
    ### Convert NB params: e.g from LJ_rmin2_eps tp LJ_sig_eps
    ###
    ConvertNonBondedParamsTo(Mol,NBFuncType)




    ###
    ### Gettting missing data (where possible) or terminate
    ###     
    D_FT={} 
    GMX_FT_Fields=['atfunc','prfunc','bdfunc','anfunc','tnfunc','imfunc']
    GMX_FT_List=['NB','NB','Bond','Angle','Torsion','Improper']
    # Must be global!
    #MolFields = ['Atoms','Pairs','Bonds','Angles','Torsions','Impropers']
    
    ### AtomNBFuncType
    lst=getattr(Mol,'Atoms')
    D_FT['Atoms']=[]
    for F in lst:
        if not F.Func == 'none':
            if not F.Func in  GMXFuncTypes['NB'].keys(): 
                sys.stderr.write('--- !Error in function '+inspect.stack()[0][3]+': the functional type ['+F.Func+'] for an Atom not present in GMXNBFuncTypes global dictionary. The molecule instance named ['+Mol.Name+']. he dictionary GMXNBFuncTypes ['+str(GMXNBFuncTypes)+']. Extend the dictionary if necessary. Aborting.\n')
                sys.exit(1)
                    
            D_FT['Atoms'].append(GMXFuncTypes['NB'][F.Func])
        else:
            sys.stderr.write('--- !Error in function '+inspect.stack()[0][3]+': in the molecule instance named ['+Mol.Name+'] the functional form for an atoms is none. The atom field is: (Aborting)\n')
            sys.stderr.write(str(F.__dict__)+'\n')
            sys.exit(1)

    # Mind! For gromacs topology all atom types must have the same functional form (it is set in [ defaults ] section)
    AtomNBFuncType = CheckIfAllElementsInArrayAreTheSameAndReturn(D_FT['Atoms'])
    if bool(AtomNBFuncType) == 0:
            sys.stderr.write('--- !Error in function '+inspect.stack()[0][3]+': in the Mol instance named ['+Mol.Name+'] the functional forms for all atoms are not identical! But this is required by Gromacs topology file format. Aborting\n')
            sys.exit(1)
        

    

    ### Set functional types for pairs
    if Vars['prfunc'] == 'none':
        # prfunc not given as the function argument. Get it from the Mol info
        
        # Fill the function types dict D_FT for PAIRS
        lst=getattr(Mol,'Pairs')
        D_FT['Pairs']=[]
        for F in lst:
                if not F.Func == 'none':
                    # If pair functional type given - use it
                    if not F.Func in  GMXFuncTypes['NB'].keys(): 
                        sys.stderr.write('--- !Error in function '+inspect.stack()[0][3]+': the functional type ['+F.Func+'] for Pair not present in GMXNBFuncTypes global dictionary. The molecule instance named ['+Mol.Name+']. he dictionary GMXNBFuncTypes ['+str(GMXNBFuncTypes)+']. Extend the dictionary if necessary. Aborting.\n')
                        sys.exit(1)
                    D_FT['Pairs'].append(GMXFuncTypes['NB'][F.Func])
                else:
                    # Set function type as in Atoms 
                    D_FT['Pairs'].append(AtomNBFuncType)
                    sys.stderr.write('--- !Warning in function '+inspect.stack()[0][3]+': in the molecule instance named ['+Mol.Name+'] the functional form for an entry in ['+'Pairs'+'] is not known (set to none). We set it to the function types of the Atom entries ['+str(AtomNBFuncType)+']. The Pair field is: (Continue)\n')
                    sys.stderr.write(str(F.__dict__))

    else:
        # prfunc not given as the function argument. Ignore the Mol info and set this to the given value.
        D_FT['Pairs'] = [ prfunc for i in range(len(Mol.Pairs))]



    
    ### Set functional types for bonded interactions
    i=1 # start from the third entry (skip Atoms and Pairs)
    for k in GMX_FT_Fields[2:]:
        i+=1
        if Vars[k] != 'none':
            # If the function types given as arguments to the function - redefine the arrays                    
            # Set the whole array to the given value
            D_FT[MolFields[i]]=[]
            D_FT[MolFields[i]] = [ Vars[k] for j in range(len(getattr(Mol,MolFields[i])))]
        else:
            # Fill the function types dict D_FT for bonded interactions
            #for k in GMX_FT_List[2:]:
                k2=GMX_FT_List[i]
                lst=getattr(Mol,MolFields[i])
                D_FT[MolFields[i]]=[]
                for F in lst:
                    if not F.Func == 'none':
                        if not F.Func in  GMXFuncTypes[k2].keys(): 
                            sys.stderr.write('--- !Error in function '+inspect.stack()[0][3]+': the functional type ['+F.Func+'] for a bonded interaction not present in GMXFuncTypes global dictionary. The molecule instance named ['+Mol.Name+']. The bobnded interaction instance ['+str(F.__dict__)+']. The dictionary GMXFuncTypes ['+str(GMXFuncTypes)+']. Extend the dictionary if necessary. Aborting.\n')
                            sys.exit(1)
                        D_FT[MolFields[i]].append(GMXFuncTypes[k2][F.Func])
                    else:
                        if Style == 'Implicit':
                            # Set default function type parameters
                            D_FT[MolFields[i]].append(GMXFuncTypes[k2]['Default'])
                            sys.stderr.write('--- !Warning in function '+inspect.stack()[0][3]+': you requested Implicit molcular topology (only atom numbers and functional types written explicitly for bonded interactions), but in the molecule instance named ['+Mol.Name+'] the functional form for an entry in ['+MolFields[i]+'] is not known (set to none). We set it to the default value ['+str(GMXFuncTypes[k2]["Default"])+'] The field is: ['+str(F.__dict__)+']. Continue\n')
                        elif Style == 'Explicit':
                            sys.stderr.write('--- !Error in function '+inspect.stack()[0][3]+': you requested Explicit molecular topology (all the parameters written explicitly), but in the molecule instance named ['+Mol.Name+'] the functional form for an entry in ['+MolFields[i]+'] is not known (set to none). The field is: ['+str(F.__dict__)+']. \n Aborting.\n')
                            sys.exit(1)
    
    ###
    ### Perform conversions
    ###

    
    # Convert units to GMX
    ConvertUnitsTo(Mol,'togmx')

    # Convert angle units to degree 
    ConvertUnitsTo(Mol,'togmx',lrad2degree='yes')




    ###
    ### Writing the ITP file
    ###
    of=open(FName,'w')

    # Header
    s=";;; This is the itp for Mol ["+Mol.Name+"] created by the ["+sys.argv[0]+"] tool "
    s2=";;; Generated by user " +user+ " at "+today
    of.write(";"*len(s)+"\n"\
            +s+"\n"\
            +s2+"\n"\
            +";"*len(s)+"\n\n")
#[ moleculetype ]
#; molname	nrexcl
#SOL		2


    # Writing FFCONV info
    string= \
    "; -- ffconv.py info: \n" 
    for k,v in Vars.iteritems():
        string += "; "+str(k)+" = "+str(v)+"\n"
    string+= "\n"
    
    # Writing Mol info
    string+= \
    "; -- molecule parameters: \n" +\
    "; FFName = "+str(Mol.FFName)+"\n" +\
    "; Units = "+str(Mol.Units)+"\n" +\
    "; GeneratedBy = "+str(Mol.GeneratedBy)+"\n" +\
    "; Comment = "+str(Mol.Comment)+"\n" +\
    "\n" 

    # Molecule name and NREXCL
    string+= \
    "[ moleculetype ]\n" +\
    "; molname     nrexcl \n" +\
    " "+Mol.Name+"    "+str(NREXCL)+"\n\n" 
    

#[ atoms ]
#;   nr   type  resnr residue  atom   cgnr     charge       mass
#     1  opls_116   1    SOL     OW      1      -0.82
#     2  opls_117   1    SOL    HW1      1       0.41
#     3  opls_117   1    SOL    HW2      1       0.41

    
    string += \
    "[ atoms ]\n"+\
    ";   nr   type  resnr residue  atom   cgnr     charge  (mass)  \n"
    for A in Mol.Atoms:
        if isnan(A.Charge): 
            sys.stderr.write('--- !Error in function '+inspect.stack()[0][3]+': the molecule instance named ['+Mol.Name+'] has atoms with charge set to NaN. The Atom instance ['+str(F.__dict__)+']. Aborting.\n ')
            sys.exit(1)
        if isnan(A.Number): 
            sys.stderr.write('--- !Error in function '+inspect.stack()[0][3]+': the molecule instance named ['+Mol.Name+'] has atoms with number set to NaN. The Atom instance ['+str(F.__dict__)+']. Aborting.\n ')
            sys.exit(1)
        if not isnan(A.Mass):
            s=" %6i %6s %6i %6s %6s %6i %14.8f %10.4f " % (A.Number, A.Type, 1, Mol.Name[0:4], A.Name, A.Number, A.Charge, A.Mass) 
        else:
            s=" %6i %6s %6i %6s %6s %6i %14.8f  " % (A.Number, A.Type, 1, Mol.Name[0:4], A.Name, A.Number, A.Charge) 
        if len(A.Quality) != 0: s += "  ; qual: "+A.Quality+""
        if len(A.Description) != 0: s += ", desc: "+A.Description+""
        s+='\n'
        string+=s
    string+='\n \n'

#[ pairs ]
#; i	j	funct	sig1-4	eps1-4
#1	2	1	0.3	0.1
#1	3	1	0.3	0.1
    Fs=Mol.Pairs
    lpnan=False
    if len(Fs) != 0:
        string += \
        "[ pairs ]\n"+\
        ";      i   j\n"
        for i,F in enumerate(Fs):
            s="%6i %6i     %6i  " % (F.AtNumbers[0],F.AtNumbers[1],D_FT['Pairs'][i]) 
            if Style == 'Explicit':
                for p in F.p:
                    if not isnan(p):
                        tmp="%16.6f" % (p)
                        s += tmp
                    else:
                        lpnan=True
                        #sys.stderr.write('--- !Warning in function '+inspect.stack()[0][3]+': you requested Explicit molecular topology, but in the molecule instance named ['+Mol.Name+'] the potential parameters for a pair 1-4 entry are NaN. Not a problem: usually one could request to generate pairs potential parameters based on the atom potential parametres if in the [ defaults ] field of the system topology you specify yes. Continue, with no output for pairs potential parameters. The pair is:\n ')
                        #sys.stderr.write(str(F.__dict__)+'\n')
                        #sys.stderr.write('--- Continue \n')
                    
            if len(F.Quality) != 0: s += "  ; qual: "+F.Quality+""
            if len(F.Description) != 0: 
                s += ", desc: "+F.Description+""
            desc = TryToGenDescription(F)
            if len(desc) !=0: s += ", "+desc+""
            s+='\n'
            string+=s

        if lpnan: 
            sys.stderr.write('--- Note in function '+inspect.stack()[0][3]+': you requested Explicit molecular topology, but in the molecule instance named ['+Mol.Name+'] SOME potential parameters for a pair 1-4 entry are NaN. Not a problem: usually one could request to generate pairs potential parameters based on the atom potential parametres if in the [ defaults ] field of the system topology you specify yes OR the missing [ pairtypes ] are given in the force field itp file explicitly. Continue, with no output for pairs with NaN potential parameters.\n ')

        string+='\n \n'

#[ bonds ]
#; i	j	funct	length	force.c.
#1	2	1	0.1	345000	0.1     345000
#1	3	1	0.1	345000	0.1     345000
    Fs=Mol.Bonds
    if len(Fs) != 0:
        string += \
        "[ bonds ]\n"+\
        ";      i   j     func     (parameters)\n"
        for i,F in enumerate(Fs):
            s="%6i %6i     %6i  " % (F.AtNumbers[0],F.AtNumbers[1],D_FT['Bonds'][i]) 
            if Style == 'Explicit':
                for p in F.p:
                    if isnan(p): 
                        sys.stderr.write('--- !Error in function '+inspect.stack()[0][3]+': you requested Explicit molecular topology, but in the molecule instance named ['+Mol.Name+'] SOME potential parameters for bonded interactions of class ['+str(F.__class__.__name__)+'] are NaN. The bonded class instance ['+str(F.__dict__)+']. Maybe try Style = Implicit ? Aborting.\n ')
                        sys.exit(1)
                    tmp="%20.6f" % (p)
                    s += tmp
            if len(F.Quality) != 0: s += "  ; qual: "+F.Quality+""
            if len(F.Description) != 0: 
                s += ", desc: "+F.Description+""
            desc = TryToGenDescription(F)
            if len(desc) !=0: s += ", "+desc+""
            s+='\n'
            string+=s
        string+='\n \n'


#[ angles ]
#; i	j	k	funct	angle	force.c.
#2	1	3	1	109.47	383	109.47	383
    Fs=Mol.Angles
    if len(Fs) != 0:
        string += \
        "[ angles ]\n"+\
        ";      i    j   k     func     (parameters)\n"

        for i,F in enumerate(Fs):
            s="%6i %6i %6i      %6i   " % (F.AtNumbers[0],F.AtNumbers[1],F.AtNumbers[2],D_FT['Angles'][i]) 
            if Style == 'Explicit':
                for p in F.p:
                    if isnan(p): 
                        sys.stderr.write('--- !Error in function '+inspect.stack()[0][3]+': you requested Explicit molecular topology, but in the molecule instance named ['+Mol.Name+'] SOME potential parameters for bonded interactions of class ['+str(F.__class__.__name__)+'] are NaN. The bonded class instance ['+str(F.__dict__)+']. Maybe try Style = Implicit ? Aborting.\n ')
                        sys.exit(1)
                    tmp="%16.6f" % (p)
                    s += tmp
            if len(F.Quality) != 0: s += "  ; qual: "+F.Quality+""
            if len(F.Description) != 0: 
                s += ", desc: "+F.Description+""
            desc = TryToGenDescription(F)
            if len(desc) !=0: s += ", "+desc+""
            s+='\n'
            string+=s
        string+='\n \n'


#[ dihedrals ]
#;  ai    aj    ak    al funct
#    9    10    12     6     2 
#    9    10    14    18     2 
    Fs=Mol.Torsions
    if len(Fs) != 0:
        string += \
        "[ dihedrals ]\n"+\
        ";      ai    aj   ak    al    func   (parameters)  \n"

        for i,F in enumerate(Fs):
            s="%6i %6i %6i %6i      %6i  " % (F.AtNumbers[0],F.AtNumbers[1],F.AtNumbers[2],F.AtNumbers[3],D_FT['Torsions'][i]) 
            if Style == 'Explicit':
                for p in F.p:
                    if isnan(p): 
                        sys.stderr.write('--- !Error in function '+inspect.stack()[0][3]+': you requested Explicit molecular topology, but in the molecule instance named ['+Mol.Name+'] SOME potential parameters for bonded interactions of class ['+str(F.__class__.__name__)+'] are NaN. The bonded class instance ['+str(F.__dict__)+']. Maybe try Style = Implicit ? Aborting.\n ')
                        sys.exit(1)
                    tmp="%16.6f" % (p)
                    s += tmp
            if len(F.Quality) != 0: s += "  ; qual: "+F.Quality+""
            if len(F.Description) != 0: 
                s += ", desc: "+F.Description+""
            desc = TryToGenDescription(F)
            if len(desc) !=0: s += ", "+desc+""
            s+='\n'
            string+=s
        string+='\n \n'

#[ dihedrals ]  ; improper
#;  ai    aj    ak    al funct
#    9    10    12     6     2 
#    9    10    14    18     2 
    Fs=Mol.Impropers
    if len(Fs) != 0:
        string += \
        "[ dihedrals ] ; improper dihedrals \n"+\
        ";      ai    aj   ak    al     func      (parameters) \n"

        for i,F in enumerate(Fs):
            s="%6i %6i %6i %6i      %6i   " % (F.AtNumbers[0],F.AtNumbers[1],F.AtNumbers[2],F.AtNumbers[3],D_FT['Impropers'][i]) 
            if Style == 'Explicit':
                for p in F.p:
                    if isnan(p): 
                        sys.stderr.write('--- !Error in function '+inspect.stack()[0][3]+': you requested Explicit molecular topology, but in the molecule instance named ['+Mol.Name+'] SOME potential parameters for bonded interactions of class ['+str(F.__class__.__name__)+'] are NaN. The bonded class instance ['+str(F.__dict__)+']. Maybe try Style = Implicit ? Aborting.\n ')
                        sys.exit(1)
                    tmp="%16.6f" % (p)
                    s += tmp
            if len(F.Quality) != 0: s += "  ; qual: "+F.Quality+""
            if len(F.Description) != 0: 
                s += ", desc: "+F.Description+""
            desc = TryToGenDescription(F)
            if len(desc) !=0: s += ", "+desc+""
            
            s+='\n'
            string+=s
        string+='\n \n'

    of.write(string)






#----------------------------------------------------------------------------------------------------------------------
def ReadGromacsMoleculeITPFile(fn):
    # TODO: complete the function 
    AtomList = [] ; AtType=[] ;  AtName=[] ; q = [] ; Mass = [] ; ResidueName=[]; ResidueNumber=[]; MolName=[]
    PairList = []
    BondList = [] ; BdNumbers = [] 
    AngleList = [] ; AnNumbers = [] 
    TorsionList = [] ; TnNumbers = [] 
    ImproperList = [] ; ImNumbers = [] 

    FieldsList=['[ moleculetype ]','[ atoms ]','[ bonds ]','[ pairs ]','[ angles ]','[ dihedrals ]']

    ifile=open(fn,'r')
    iBond = 0 ; iAngle = 0 ; iTorsion = 0 ; iImproper=0
    while 1:
        line=ifile.readline()
    
        if not line: break
        ln=line.split() ; 
        if len(ln)<1: continue
        ln=line.split(';'); ln=ln[0] # skip comments
        
        for val in FieldsList:
            field='none'
            if (ln.find(val)!=-1):
                field=ln[1]
                break

#[ moleculetype ]
#; name  nrexcl
#NMP   3
        if   field=='[ moleculetype ]':
            while ():
                line=ifile.readline()
                if not line: break
                ln=line.split() ; 
                if len(ln)<1: continue
                ln=line.split(';'); ln=ln[0] # skip comments
                ln=line.split() ; 
                if (len(ln)<1): continue
                if (ln.find('[') != -1 ):
                    ifile.seek(ifile.tell() - len(l))
                    break
                MoleculeName=ln[0]
                nrexcl=int(ln[2])

#[ atoms ]
#;   nr  type    resnr   residu  atom    cgnr    charge  mass
#   1   C1    1     NMP   C11    1    -0.05000    12.011  ; -0.110    0.0261950000   
#   2   C1    1     NMP   C12    2    -0.12000    12.011  ; -0.110   -0.0512900000   
#   3   C1    1     NMP   C13    3    -0.12000    12.011  ; -0.110   -0.0782780000   
        elif field=='[ atoms ]':
            while ():
                line=ifile.readline()
                if not line: break
                ln=line.split() ; 
                if len(ln)<1: continue
                ln=line.split(';'); ln=ln[0] # skip comments
                ln=line.split() ; 
                if len(ln)<1: continue
                if (ln.find('[') != -1 ):
                    ifile.seek(ifile.tell() - len(l))
                    break
                AtomList.append(int(ln[0]))
                AtType.append(ln[1])
                ResidueNumber.append(int(ln[2]))
                ResidueName.append(ln[3])
                AtName.append(ln[4])
                # skipping the charge group. this is old fashioned and is not used in the last force fields
                q.append(float(ln[6]))
                Mass.append(float(ln[7]))

### From this point this is not finished.
#[ pairs ]
#;   ffio_pairs[70] { i_ffio_ai i_ffio_aj s_ffio_funct r_ffio_c1 :::
#      1     8     1   
#      1     16    1   
        elif field=='[ bonds ]':
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

        ### Not finished! Old code here !

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
    prNumbers=[]
    Mol.FillMoleculeTopology(AtType,AtName,q,BdNumbers,PrNumbers,AnNumbers,TnNumbers,ImNumbers)
    
    #print Mol
    return Mol



def TryToGenDescription(O):
    desc=''
    if hasattr(O,'AtNames'):
        if not 'none' in getattr(O,'AtNames'):
            desc += str(getattr(O,'AtNames'))
        
    if hasattr(O,'AtTypes'):
        if not 'none' in getattr(O,'AtTypes'):
            s = str(getattr(O,'AtTypes'))
            if len(desc) != 0: 
                desc += ', '+s
            else:
                desc += s
    return desc
        
    

  
    
    



