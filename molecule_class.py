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

# Last Modified: May 2014

from copy import copy
import math
import inspect
import xml.etree.ElementTree
import traceback
from FactorToConvertUnits import *
from elements import *

from datetime import *
import getpass

# TODO
# - Documentation for each function
# - put sys.stderr.write to all error / warning messages
# - add traceback to all Error messages



###
### Set global variables
### 
user = str(getpass.getuser())
today = datetime.today()
NaN = float('NaN')


# Global variables
FF_OPLS_family = ['OPLSAA','OPLSUA','OPLSAA_NONBOND_ONLY']
FF_CHARMM_family = ['CHARMM27','CGENFF','CGENFF2B7','CGENFF2B7_NONBOND_ONLY']
FF_AMBER_family = ['GAFF','GAFF_NONBOND_ONLY']

# Defaults for FFs:
# - fudgeLJ - factor to sclae LJ parameters for 1-4 interactions
# - fudgeQQ - factor to sclae Coulomb interaction for 1-4 interactions
# - NREXCL - number of nearest bonded atom to put exclusions list
# - CombRule - combining rules 
FFvars = {'CHARMM': 
                { 'FFNames': FF_CHARMM_family,                
                  'fudgeLJ':1.0 , 'fudgeQQ':1.0, 'NREXCL':3, 'CombRule':'LorentzBerthelot' },
          'OPLS':
                { 'FFNames': FF_OPLS_family,                
                  'fudgeLJ':0.5 , 'fudgeQQ':0.5, 'NREXCL':3, 'CombRule':'Geometric' },
          'AMBER':
                { 'FFNames': FF_AMBER_family,                
                  'fudgeLJ':0.5 , 'fudgeQQ':0.833333, 'NREXCL':3, 'CombRule':'LorentzBerthelot' }
}

KnownFFs = []
for FFclass,val in FFvars.iteritems():
    for val2 in val['FFNames']:
        KnownFFs.append(val2) 
KnownFFFamilies = FFvars.keys()

def GetFFFamily(FFName):
    for FFFamily,val in FFvars.iteritems():
        if FFName.upper() in val['FFNames']:
            return FFFamily
    sys.stderr.write("--- !Warning in "+inspect.stack()[0][3]+": given FFName ["+FFName+"] is not present in FFvars dictionary (meaning it is not known). Maybe extend the FFvars ? Now FFvars: ["+str(FFvars)+"]. Returning none. \n")
    return 'none'
    

KnownAtFunc = ['LJ_sig_eps','LJ_Rmin2_eps','LJ_A_B','Buckingham']
KnownClassNames = ['Atom','Pair','Bond','Angle','Torsion','Improper','AtomType','BondType','PairType','AngleType','TorsionType','ImproperType']
KnownAttributes=['AtNames','AtNumbers','List','AtTypes','Func','p','p_14','Name','Type','Mass','PeriodicTableNum','X','Y','Z','Charge']

# Here the fields in the lists correspond to each other. The order matters.
MolFields = ['Atoms','Pairs','Bonds','Angles','Torsions','Impropers']
FF_Fields = ['AtomTypes','PairTypes','BondTypes','AngleTypes','TorsionTypes','ImproperTypes']
FT_Fields = ['NB','NB','Bond','Angle','Torsion','Improper']





###--------------------------------------------------------------------------------------------------------------------------------
### Some auxiliary functions and general functions 

def UnitsForFuncTypes(Units):
    
    E=Units['Energy']
    D=Units['Distance']
    A=Units['Angle']
    N=Units['Number']

    EN=E+'/'+N
    Dsq=D+'^2.0'
    Asq=A+'^2.0'

    NBFuncTypes = {'LJ_sig_eps': [D,EN], 'LJ_Rmin2_eps': [D,EN], 'Buckingham':[]}
    BondFuncTypes = {'Harmonic': [D,EN+'/'+Dsq]} # To be extended
    AngleFuncTypes = {'Harmonic': [A,EN+'/'+Asq], 'UreyBradley': [A,EN+'/'+Asq,D,EN+'/'+Dsq]} # To be extended
    TorsionFuncTypes = {'Periodic': [A,EN,1], 'PeriodicMultiple':[A,EN,1], 'Fourier':[EN,EN,EN,EN,EN,EN], 'RB':[EN,EN,EN,EN,EN,EN]} # To be extended
    ImproperFuncTypes = {'Periodic': [A,EN,1], 'PeriodicMultiple':[A,EN,1], 'Harmonic': [A,EN+'/'+Asq], 'Fourier':[EN,EN,EN,EN,EN,EN], 'RB':[EN,EN,EN,EN,EN,EN] } # To be extended

    FuncTypes={}
    FuncTypes['NB']=NBFuncTypes
    FuncTypes['Bond']=BondFuncTypes
    FuncTypes['Angle']=AngleFuncTypes
    FuncTypes['Torsion']=TorsionFuncTypes
    FuncTypes['Improper']=ImproperFuncTypes
    
    #print NBFuncTypes.keys() 
    #print FuncTypes
    #raw_input()

    return FuncTypes





def ApplyCombRules(p1,p2,CombRule,Func):

    KnownFT=['LJ_sig_eps','LJ_Rmin2_eps','LJ_A_B']
    if not Func in KnownFT:
        sys.stderr.write("--- !Error in "+inspect.stack()[0][3]+": nonbonded unctional type ["+str(Func)+"] not known. Cannot calculate mixed potential parameters. iKnown funct types: ["+str(KnownFT)+"]. Aborting.\n")
        sys.exit(1)
        
    KnownCR=['LorentzBerthelot','Geometric']
    if not CombRule in KnownCR:
        sys.stderr.write("--- !Error in "+inspect.stack()[0][3]+": combining rule not known ["+str(CombRule)+"]. Cannot estimate mixed 1-4 nonbonded parameters. Known rules: ["+str(KnownCR)+"]. Aborting.\n")
        sys.exit(1)
    
    if Func == 'LJ_A_B':
        # Tricky potential type. Convert functiontype first to LJ_sig_eps, evaluate mixed parameters, convert mixed parameter to inital func type 
        p1 = ConvertAtomFuncParam(p1,Func,'LJ_sig_eps')        
        p2 = ConvertAtomFuncParam(p1,Func,'LJ_sig_eps')
                
    s1=p1[0];e1=p1[1]
    s2=p2[0];e2=p2[1]
    if CombRule == 'LorentzBerthelot':
        s12 = (s1+s2)/2.0 ; e12 = math.sqrt(e1*e2)
    elif CombRule == 'Geometric':
        s12 = math.sqrt(s1*s2) ; e12 = math.sqrt(e1*e2)
        
    if Func == 'LJ_A_B':    
        [s12,e12] = ConvertAtomFuncParam([s12,e12],'LJ_sig_eps',Func)
        
    return [s12,e12]


        
def FindObjByFieldValue(ObjList,Fname,Fval):
            toreturn=[]
            for obj in ObjList:
                if Fval == getattr(obj,Fname):
                    toreturn.append(copy(obj))
            
            if len(toreturn) < 1:
                sys.stderr.write("--- !Error in "+inspect.stack()[0][3]+": there is no entries in the object list where attribute ["+Fname+"] has value ["+str(Fval)+"].\n")
                sys.exit(1)
            elif len(toreturn) > 1:
                sys.stderr.write("--- !Warning in "+inspect.stack()[0][3]+": there are multiple entries in the object list where attribute ["+Fname+"] has value ["+str(Fval)+"]. Returning the first instance.\n")
            toreturn=toreturn[0]
            return toreturn



def cpAllAttributes(O1,O2,excl_list):
# This function is probably redundant and can be removed by python standart functions. Check documentation.
    for attr, value in O1.__dict__.iteritems():
        if not attr in excl_list:
            setattr( O2, attr, copy(value))
    return



def FillMasses(O): 
    if   hasattr(O,'Atoms'):
        l=getattr(O,'Atoms')
    elif hasattr(O,'AtomTypes'):    
        l=getattr(O,'AtomTypes')
    else:
        sys.stderr.write("--- !Warning in "+inspect.stack()[0][3]+": the input object has no attributes Atoms or Atomypes. Nothing to do. Continue.\n")
        return
    
    for i,v in enumerate(l):
        try:
            v.Mass = copy(Elements[TrimElementName(v.Type)]['Mass'])
            v.PeriodicTableNum = copy(Elements[TrimElementName(v.Type)]['PeriodicTableNum'])
        except:
            sys.stderr.write("--- !Warning in "+inspect.stack()[0][3]+": failed to identify element for string ["+str(v.Type)+"] object ["+str(v.__dict__)+"]. Maybe extend the Elements dictionary [see elements.py] ?. Continue.\n")
    return    



def ConvertAtomFuncParam(p1,Func1,Func2):
    p2=[NaN,NaN] # initialize
    # Convert first to 'LJ_sig_eps', then from this to a different form
    if Func1 == 'LJ_sig_eps':
        # No conversion needed
        p2=p1
    elif Func1 == 'LJ_Rmin2_eps':
        p2[0] = p1[0] * 2.0 / 2.0**0.1666666666666
        p2[1] = p1[1]
    elif Func1 == 'LJ_A_B':
        A=p1[0] ; B=p1[1]
        p2[0] = ( A / B )**0.166666666666666
        p2[1] = B**2 / 4 / A 
    else:
        if Func1 in KnownAtFunc:
            sys.stderr.write("--- !Error in "+inspect.stack()[0][3]+": function type is known, but the function doesnot know how to convert this to LJ_sig_eps yet. The code has to be extended. Aborting.\n")
            sys.exit(1)
        else:
            sys.stderr.write("--- !Warning in "+inspect.stack()[0][3]+": function type ["+Func1+"] NOT known. Check for mistakes, or develop the code. Aborting.\n")
            sys.exit(1)
    
    # Convert from LJ_sig_eps to Func2
    if Func2 == 'LJ_sig_eps':
        p2=p2
    elif Func2 == 'LJ_Rmin2_eps':
        p2[0] = p2[0] / ( 2.0 / 2.0**0.1666666666666 )
        p2[1] = p2[1]
    elif Func2 == 'LJ_A_B':
        sig=p2[0]; eps=p2[1]
        p2[0] = 4.0 * eps * sig**12.0
        p2[1] = 4.0 * eps * sig**6.0
    else:
        sys.stderr.write("--- !Error in "+inspect.stack()[0][3]+": do not know how to convert from internal form LJ_sig_eps to this functional form ["+str(Func2)+"]. Extend the code. Aborting.\n")
        sys.exit(1)
    
    return p2




def ConvertNonBondedParamsTo(O,toFunc):
    # Convert first to 'LJ_sig_eps'
    
    if   O.__class__.__name__ == 'ForceField':
        for F in O.AtomTypes:
            setattr( F, 'p', ConvertAtomFuncParam(F.p, F.Func, toFunc)  ) 
            setattr( F, 'p_14', ConvertAtomFuncParam(F.p_14, F.Func, toFunc)  ) 
            setattr( F, 'Func', toFunc ) 
        for F in O.PairTypes:
            setattr( F, 'p', ConvertAtomFuncParam(F.p, F.Func, toFunc)  ) 
            setattr( F, 'Func', toFunc ) 
        
    elif O.__class__.__name__ == 'Molecule':
        for F in O.Atoms:
            setattr( F, 'p', ConvertAtomFuncParam(F.p, F.Func, toFunc)  ) 
            setattr( F, 'p_14', ConvertAtomFuncParam(F.p_14, F.Func, toFunc)  ) 
            setattr( F, 'Func', toFunc ) 
        for F in O.Pairs:
            setattr( F, 'p', ConvertAtomFuncParam(F.p, F.Func, toFunc)  ) 
            setattr( F, 'Func', toFunc ) 
        
    else:
        print 'Error !!!: function '+inspect.stack()[0][3]+' doesnot know how to work with class '+O.__class__.__name__+'. Check for mistakes or extend the code. Exit.' 
        sys.exit(1)
    return





def ConvertUnitsTo(obj,unitkey,**kwargs):

    Vars={}
    Vars['lrad2degree']='no'
    
    ###
    ### Process function arguments 
    ###
    for k,v in kwargs.iteritems():
        if not k in Vars.keys():
            sys.stderr.write('--- !Warning in function '+inspect.stack()[0][3]+' not known input key argument '+k+'. Continue.\n')
        for key,val in Vars.iteritems():
            if k == key: Vars[k]=v

    ###
    ### Checks 
    ###
    if Vars['lrad2degree'] == 'yes':
        if obj.Units != SpecialKeyUnits[unitkey]:
            sys.stderr.write('--- !Warning in function '+inspect.stack()[0][3]+' you want to convert all angles from rad to degree, but object units ['+str(obj.Units)+'] and unit to convert to ['+str(SpecialKeyUnits[unitkey])+'] are NOT identical. To run the function with the lrad2degree=yes key one first should convert units without lrad2degree=yes key. And only after, convert angles from rad to degree. No converion made. Return.\n')
            return

    FuncTypes = UnitsForFuncTypes(obj.Units)

    def CalcFuncTypesFactors(FT):
       # FT - funct types
        FTfactors={}
        for k,v in FT.iteritems():
            FTfactors[k]={}
            for k2,v2 in v.iteritems():
                tmp=[]
                for v3 in v2:
                    if not isinstance(v3,int):
                        if Vars['lrad2degree'] == 'yes' and v3 == 'rad':
                            tmp.append(FactorToConvertUnits(v3,'degree'))
                        else:
                            tmp.append(FactorToConvertUnits(v3,unitkey))
                            
                    else:
                        # E.g. multiplicity - no need to convert
                        tmp.append(v3)
                FTfactors[k][k2]=tmp
        return FTfactors                    
    


    FTfactors = CalcFuncTypesFactors(FuncTypes) 
    

    if obj.__class__.__name__ == 'ForceField':
        ObjFields = FF_Fields
    elif obj.__class__.__name__ == 'Molecule':
        ObjFields = MolFields
    else:
        # TODO Not known
        pass
        
    for i,FN in enumerate(ObjFields):
        Fs=getattr(obj,FN) 
        for F in Fs:
            if F.Func == 'none':
                continue
            if not F.Func in FuncTypes[FT_Fields[i]].keys():
                sys.stderr.write("--- !Warning in "+inspect.stack()[0][3]+": an object functional type ["+str(F.Func)+"] is not known AND not set to none. The object class: ["+str(F.__class__.__name__)+"], object: ["+str(F.__dict__)+"]. Maybe extend the FTfactors dictionary?. Now it is ["+str(FTfactors)+"] Continue.\n")
                
            c = FTfactors[FT_Fields[i]][F.Func]
            # Mind! If there are no parameters given len(F.p) can be != len(c). But len(F.p) always smaller
            for j in range(len(F.p)):
                F.p[j] *= c[j]
                try:    
                    F.p_14[j] *= c[j]
                except:
                    pass
    
    if Vars['lrad2degree'] != 'yes':
        obj.Units = copy(SpecialKeyUnits[unitkey.lower()])
    else:
        obj.Comment += ' !!! Mind all angles were converted from RAD to DEGREES, but angles in the force constants are in RAD!'

    return
            



#======================================================================================================================
#======================================================================================================================
###--------------------------------------------------------------------------------------------------------------------
### Description of classes

class AtomType:
    def __init__(self):
        self.Type = 'none'
        self.Mass = NaN
        self.PeriodicTableNum = 0
        #self.Charge = 'NaN'
        #self.Func = 'LJ_sig_eps' # LJ_Rmin2_eps, LJ_sig_eps, 
        self.Func = 'none' # LJ_Rmin2_eps, LJ_sig_eps, 
        self.p = [NaN,NaN]
        self.p_14 = [NaN,NaN] # for [ pairs ], 1-4 pair interactions. 
        self.Description = ''
class Atom(AtomType):
    def __init__(self,AtNum):
        AtomType.__init__(self)
        #self.__data__ = AtTy
        self.Number = AtNum
        self.X = NaN ; self.Y = NaN ; self.Z = NaN
        self.Name = 'none'
        self.Charge = NaN
        self.Quality = 'none'
        
class BondType:
    def __init__(self):
        self.AtTypes = ['none','none']
        #self.Func = 'Harmonic' # Harmonic 
        self.Func = 'none' # Harmonic 
        self.p = [NaN,NaN]
        self.Description = ''
class Bond(BondType):
    def __init__(self,BdNum):
        BondType.__init__(self)
        self.Number = BdNum
        self.AtNames = ['none','none']
        self.AtNumbers = [NaN,NaN]
        self.Quality = 'none'

        
class PairType:
    def __init__(self):
        self.AtTypes = ['none','none']
        #self.Func = 'LJ_sig_eps' # Harmonic 
        self.Func = 'none' # Harmonic 
        self.p = [NaN,NaN]
        self.Description = ''
class Pair(PairType):
    def __init__(self,PrNum):
        PairType.__init__(self)
        self.Number = PrNum
        self.AtNames = ['none','none']
        self.AtNumbers = [NaN,NaN]
        self.Quality = 'none'



class AngleType:
    def __init__(self):
        self.AtTypes = ['none','none','none']
        #self.Func = 'Harmonic' # Harmonic, UreyBradley
        self.Func = 'none' # Harmonic, UreyBradley
        self.p = [NaN,NaN]
        self.Description = ''
class Angle(AngleType):
    def __init__(self,AnNum):
        AngleType.__init__(self)
        self.Number = AnNum
        self.AtNames = ['none','none','none']
        self.AtNumbers = [0,0,0]
        self.Quality = 'none'


class TorsionType:
    def __init__(self):
        self.AtTypes = ['none','none','none','none']
        #self.Func = 'Periodic' # Periodic, PeriodicMultiple
        self.Func = 'none' # Periodic, PeriodicMultiple
        self.p = [NaN,NaN,NaN]
        self.Description = ''
        self.Quality = 'none'
        
class Torsion(TorsionType):
    def __init__(self,TnNum):
        TorsionType.__init__(self)
        self.Number = TnNum
        self.AtNames = ['none','none','none','none']
        self.AtNumbers = [0,0,0,0]
        self.Quality = 'none'


class ImproperType:
    def __init__(self):
        self.AtTypes = ['none','none','none','none']
        #self.Func = 'Harmonic' # Harmonic
        self.Func = 'none' # Harmonic
        self.p = [NaN,NaN]
        self.Description = ''
        
class Improper(ImproperType):
    def __init__(self,TnNum):
        ImproperType.__init__(self)
        self.Number = TnNum
        self.AtNames = ['none','none','none','none']
        self.AtNumbers = [0,0,0,0]
        self.Quality = 'none'

class Exclusion:
    def __init__(self,Num):
        self.Number = Num
        self.AtNames = ['none','none']
        self.AtNumbers = [0,0]
        self.Description = ''



#======================================================================================================================
#======================================================================================================================
class ForceField:
    
    def __init__(self,Name):
        self.Units = {}
        self.AtomTypes = []
        self.BondTypes = []
        self.AngleTypes = []
        self.TorsionTypes = []
        self.PairTypes = []
        self.ImproperTypes = []
        self.Name = Name.upper()
        self.Comment=''
        self.GeneratedBy=''

        for FFclass,val in FFvars.iteritems():
            if Name.upper() in val['FFNames']:
                self.fudgeLJ = val['fudgeLJ']
                self.fudgeQQ = val['fudgeQQ']
                self.NREXCL = val['NREXCL']
                self.CombRule = val['CombRule']
                self.Family = FFclass

        if not Name.upper() in KnownFFs:
            sys.stderr.write("--- !Warning in "+inspect.stack()[0][3]+": the FF name ["+Name+"] is not known. Set the fudgeLJ, fudgeQQ, NREXCL to NaN and ComRule to 'none'. Give another name or extend the list of known FFs: "+str(KnownFFs)+". Continue.\n")
            self.fudgeLJ = NaN
            self.fudgeQQ = NaN
            self.NREXCL = NaN
            self.CombRule = 'none'
            self.Family = 'none'
            

    def __repr__(self):

        # Writing FF info
        string= \
        " -- force field parameters: \n" +\
        " Name = "+str(self.Name)+"\n" +\
        " Family = "+str(self.Family)+"\n" +\
        " Units = "+str(self.Units)+"\n" +\
        " GeneratedBy = "+str(self.GeneratedBy)+"\n" +\
        " Comment = "+str(self.Comment)+"\n" +\
        " NREXCL = "+str(self.NREXCL)+"\n" +\
        " fudgeLJ = "+str(self.fudgeLJ)+"\n" +\
        " fudgeQQ = "+str(self.fudgeQQ)+"\n" +\
        " CombRule = "+str(self.CombRule)+"\n" +\
        "\n" 
        print string
        
        #from pprint import pprint
        print "-"
        print "ATOM TYPES"
        for i in range(len(self.AtomTypes)):
            print "i:"+str(i), self.AtomTypes[i].__dict__
        print "-"
        print "PAIR 14 TYPES"
        for i in range(len(self.PairTypes)):
            print "i:"+str(i), self.PairTypes[i].__dict__
        print "-"
        print "BOND TYPES"
        for i in range(len(self.BondTypes)):
            print "i:"+str(i), self.BondTypes[i].__dict__
        print "-"
        print "ANGLE TYPES"
        for i in range(len(self.AngleTypes)):
            print "i:"+str(i), self.AngleTypes[i].__dict__
        print "-"
        print "TORSION TYPES"
        for i in range(len(self.TorsionTypes)):
            print "i:"+str(i), self.TorsionTypes[i].__dict__
            #pprint(self.TorsionTypes[i].__dict__,indent=4)
        print "-"
        print "IMPROPER TYPES"
        for i in range(len(self.ImproperTypes)):
            print "i:"+str(i), self.ImproperTypes[i].__dict__
            #pprint(self.TorsionTypes[i].__dict__,indent=4)
        return "--------" 


    def GenPairParamFromAtomTypeParam(self):
        import math
        if self.CombRule == 'none' or self.CombRule == '':
            sys.stderr.write("--- !Error in "+inspect.stack()[0][3]+": the FF conbining rules must be specified. Now: ["+str(self.CombRule)+"]\n")
            sys.exit(1)
        
        
        #pairs14=[[float('NaN')]]*
        pairs14=[]
        N = len(self.AtomTypes)
        for i in range(N):
            for j in range(i,N):
                s1=self.AtomTypes[i].p_14[0] ; e1=self.AtomTypes[i].p_14[1]
                s2=self.AtomTypes[j].p_14[0] ; e2=self.AtomTypes[j].p_14[1]
                Func = CheckIfAllElementsInArrayAreTheSameAndReturn([self.AtomTypes[i].Func,self.AtomTypes[j].Func])
                if not (isnan(s1) or isnan(e1) or isnan(s2) or isnan(e2)):
                    s12,e12 = ApplyCombRules([s1,e1],[s2,e2],self.CombRule,Func)
                else:
                    sys.stderr.write("--- !Warning in "+inspect.stack()[0][3]+": some AtomType.p_14 parameters are NaN for atomtypes: ["+str(self.AtomTypes[i].Type)+"] and ["+str(self.AtomTypes[j].Type)+"]. Set the pair parameters to NaN as well. Continue.\n")
                    s12=NaN ; e12=NaN
                pairs14.append([s12,e12])
        
        tmp=[]
        cnt=-1
        for i in range(N):
            for j in range(i,N):
                cnt+=1
                AT1=self.AtomTypes[i] ; AT2=self.AtomTypes[j]
                try:
                    tmp.append(self.PairTypes[cnt])
                    sys.stderr.write("--- !Warning in "+inspect.stack()[0][3]+": PairType numbered with "+str(cnt)+" exists. It will be overwritten by the function: "+inspect.stack()[0][3]+". The overwritten PairTypes returned as the second output.\n")
                except:
                    self.PairTypes.append( PairType() )
                
                self.PairTypes[cnt].AtTypes = [ AT1.Type , AT2.Type ]
                self.PairTypes[cnt].p = copy(pairs14[cnt])
                self.PairTypes[cnt].Func = AT1.Func
                self.PairTypes[cnt].Description = 'ffconv generated'
         
        return tmp   





        
    def RemoveIdenticalFFEntries(self):
        def RemoveIdenticalFFEntriesFromList(typelist,field):
            #from operator import isSequenceType
            import types
            uniquelist=[]
            uniquelist_att=[]
            for entry1 in typelist:
                    l1=getattr(entry1,field)
                    
                    if isinstance(l1, types.ListType):
                        if not l1 in uniquelist_att and not l1[::-1] in uniquelist_att:
                            # entry is unique
                            uniquelist.append(copy(entry1))
                            uniquelist_att.append(l1)
                    else:
                        if not l1 in uniquelist_att:
                            # entry is unique
                            uniquelist.append(copy(entry1))
                            uniquelist_att.append(l1)

            return uniquelist
        
        self.AtomTypes = copy(RemoveIdenticalFFEntriesFromList(self.AtomTypes,'Type'))
        self.PairTypes = copy(RemoveIdenticalFFEntriesFromList(self.PairTypes,'AtTypes'))
        self.BondTypes = copy(RemoveIdenticalFFEntriesFromList(self.BondTypes,'AtTypes'))
        self.AngleTypes = copy(RemoveIdenticalFFEntriesFromList(self.AngleTypes,'AtTypes'))
        self.TorsionTypes = copy(RemoveIdenticalFFEntriesFromList(self.TorsionTypes,'AtTypes'))
        self.ImproperTypes = copy(RemoveIdenticalFFEntriesFromList(self.ImproperTypes,'AtTypes'))
           
    def addFF(self,FFextra):
       # TODO: checks, remove identical ? etc. 
        self.AtomTypes += FFextra.AtomTypes
        self.PairTypes += FFextra.PairTypes
        self.BondTypes += FFextra.BondTypes
        self.AngleTypes += FFextra.AngleTypes
        self.TorsionTypes += FFextra.TorsionTypes
        self.ImproperTypes += FFextra.ImproperTypes


    def CheckIfTheSameAtFuncAndReturn(self):
        i=0
        for AtomType in self.AtomTypes:
            if i==0:
                func0=AtomType.Func
            else:
                if func0 != AtomType.Func:
                    print "Warning AtomType.Func are not the same for all atoms !!!"
                    return 'none'
            i+=1
            return func0
    
    def ChangeJointAtomTypeToAtomTypesBelongingToIt(self):
    ### The function changes atomtypes with '*' to the atom types belonging tot, eg 'CA*' to 'CA1, CA2, ...'
	from copy import deepcopy
        
	ListAllAtomTypes = []
        for val in self.BondTypes + self.AngleTypes + self.TorsionTypes + self.ImproperTypes:
            for AT in val.AtTypes:
                if AT not in ListAllAtomTypes: ListAllAtomTypes.append(AT)
        
        InitFFAtomTypes = self.AtomTypes
        self.AtomTypes= []

        j=-1
        for val in InitFFAtomTypes:
	    AT=val.Type
            tmp=j 
	    for AT2 in ListAllAtomTypes:
	        if AT2 == AT:
                    j+=1
                    # Just add the AT
	            self.AtomTypes.append( copy(val) )
	        elif AT2[:-1] == AT[:-1] and AT[-1] == '*':
                    j+=1
                    # Add
		    self.AtomTypes.append( copy(val) )
		    self.AtomTypes[-1].Type= AT2
	    if tmp == j:	    
	        j+=1
	        #print '!!! Warning : atom type in the FF file is not used in the FF bonds angles, etc.'
	        self.AtomTypes.append( copy(val) )





    def FindAtomTypeHavingStringAndReturn(self,ATstr):
        cnt=0
        for AT in self.AtomTypes:
            if AT.Type == ATstr:
                #print AT.Type, ATstr
                ATreturn = AT
                cnt+=1
        if cnt < 1:
            print '!!! Error: no specified atomtypes is found in FF. Input AT: '+ATstr
            return 1
        elif cnt > 1:
            print '!!! Warning: there are more than 1 atomtypes in FF with the same name: '+Atstr+'. The last by number is taken.'
            return ATreturn
        else:
            return ATreturn



#======================================================================================================================
#======================================================================================================================
class Molecule:
    def __init__(self,Name):
        self.Units = {}
        self.Atoms = []
        self.Bonds = []
        self.Pairs = []
        self.Angles = []
        self.Torsions = []
        self.Impropers = []
        self.Exclusions = []
        self.Name=Name
        self.FFName='none'
        self.GeneratedBy='none'
        self.Comment=''
        self.ExtraData={}
        
    def __repr__(self):
    # Not changed from FF class    
        # Writing Molecule info
        string= \
        " -- force field parameters: \n" +\
        " Units = "+str(self.Units)+"\n" +\
        " FFName = "+str(self.FFName)+"\n" +\
        " Name = "+str(self.Name)+"\n" +\
        " GeneratedBy = "+str(self.GeneratedBy)+"\n" +\
        " Comment = "+str(self.Comment)+"\n" +\
        " ExtraData = "+str(self.ExtraData)+"\n" +\
        "\n" 
        print string

        from pprint import pprint
        print "-"
        print "ATOMS"
        for i in range(len(self.Atoms)):
            print "i:"+str(i), self.Atoms[i].__dict__
        print "-"
        print "PAIRS14"
        for i in range(len(self.Pairs)):
            print "i:"+str(i), self.Pairs[i].__dict__
        print "-"
        print "BONDS"
        for i in range(len(self.Bonds)):
            print "i:"+str(i), self.Bonds[i].__dict__
        print "-"
        print "ANGLES"
        for i in range(len(self.Angles)):
            print "i:"+str(i), self.Angles[i].__dict__
        print "-"
        print "TORSIONS"
        for i in range(len(self.Torsions)):
            print "i:"+str(i), self.Torsions[i].__dict__
        print "-"
        print "IMPROPERS"
        for i in range(len(self.Impropers)):
            print "i:"+str(i), self.Impropers[i].__dict__
        print "EXCLUSIONS"
        for i in range(len(self.Exclusions)):
            print "i:"+str(i), self.Exclusions[i].__dict__
        return "--------" 

    def FillAtomsFF(self,FF):
        N = len(self.Atoms)
        for i in range(N):
            #FF.AtomTypes
            AtNum=i
            AT = FF.FindAtomTypeHavingStringAndReturn(self.Atoms[AtNum].Type)            
            #self.Atoms[AtNum].cpAllAtomTypeAttributes(AT)
            cpAllAttributes(AT , self.Atoms[AtNum], ['Charge'] )
    

    def GenConnectMap(self):
        connect_map={}
        for i in range(len(self.Atoms)):
            connect_map[i+1]=[]
           
        for BD in self.Bonds:
            connect_map[BD.AtNumbers[0]].append(BD.AtNumbers[1])
            connect_map[BD.AtNumbers[1]].append(BD.AtNumbers[0])
        return connect_map
        

    def GenPairsFromBondConnectivity(self):
        pairs14=[]
        
        connect_map = self.GenConnectMap()

        for BD in self.Bonds:
            for i in connect_map[BD.AtNumbers[0]]:
                for j in connect_map[BD.AtNumbers[1]]:
                    if i != BD.AtNumbers[1] and j != BD.AtNumbers[0]: pairs14.append([i,j])
        
        pairs14_tofix = pairs14
        pairs14 = []
        for pair in pairs14_tofix:
            if not bool(set(connect_map[pair[0]]) & set(connect_map[pair[1]])): pairs14.append(pair)
        

        # Check if pairs are unique. The ring structures (e.g. benzene ring) result in double account of some 1-4 pairs. We need to exclude them.

        unique_pairs14 = CheckIfUniqueListAndReturn(pairs14)
        
        ### Pairs
        if self.Pairs != []:
            self.Pairs = []
            sys.stderr.write("--- !Warning in "+inspect.stack()[0][3]+": Mol.Pairs are not empty! Erase them and fill new ones according to bond connectivity. Continue.\n")
        for i,Pr in enumerate(unique_pairs14):
            P = Pair(i+1)
            P.AtNumbers = copy(Pr)
            P.Description = 'generated by ffconv'
            
            # get the pair functional type based on the functional type of atoms 
            A1 = FindObjByFieldValue(self.Atoms,'Number',P.AtNumbers[0])
            A2 = FindObjByFieldValue(self.Atoms,'Number',P.AtNumbers[1])
            P.Func = CheckIfAllElementsInArrayAreTheSameAndReturn( [ A1.Func , A2.Func ])
            
            self.Pairs.append(copy(P))

                    


    def GenAnglesTorsionsFromBondConnectivity(self):
        torsions=[]
        torsions_AN=[]
        torsions_AT=[]
        angles=[]        
        angles_AN=[]        
        angles_AT=[]        
        connect_map = self.GenConnectMap()
        # TODO rewrite with FindObj
        for BD in self.Bonds:
            for i in connect_map[BD.AtNumbers[0]]:
                for j in connect_map[BD.AtNumbers[1]]:
                    if i != BD.AtNumbers[1] and j != BD.AtNumbers[0]:
                        A1 = FindObjByFieldValue(self.Atoms, 'Number', i)
                        A2 = FindObjByFieldValue(self.Atoms, 'Number', BD.AtNumbers[0])
                        A3 = FindObjByFieldValue(self.Atoms, 'Number', BD.AtNumbers[1])
                        A4 = FindObjByFieldValue(self.Atoms, 'Number', j)
                        #torsions.append([i,BD.AtNumbers[0],BD.AtNumbers[1],j])
                        #torsions_AN.append([self.Atoms[i-1].Name,self.Atoms[BD.AtNumbers[0]-1].Name,self.Atoms[BD.AtNumbers[1]-1].Name,self.Atoms[j-1].Name])
                        #torsions_AT.append([self.Atoms[i-1].Type,self.Atoms[BD.AtNumbers[0]-1].Type,self.Atoms[BD.AtNumbers[1]-1].Type,self.Atoms[j].Type])
                        torsions.append([i,BD.AtNumbers[0],BD.AtNumbers[1],j])
                        torsions_AN.append([A1.Name,A2.Name,A3.Name,A4.Name])
                        torsions_AT.append([A1.Type,A2.Type,A3.Type,A4.Type])

        for key,val in connect_map.iteritems():
            for i,v1 in enumerate(val):
                for j,v2 in enumerate(val[i+1:]):
                    A1 = FindObjByFieldValue(self.Atoms, 'Number', v1)
                    A2 = FindObjByFieldValue(self.Atoms, 'Number', key)
                    A3 = FindObjByFieldValue(self.Atoms, 'Number', v2)
                    angles.append([v1,key,v2])
                    angles_AN.append([A1.Name,A2.Name,A3.Name])
                    angles_AT.append([A1.Type,A2.Type,A3.Type])
                    #angles_AN.append([self.Atoms[v1-1].Name,self.Atoms[key-1].Name,self.Atoms[v2-1].Name])
                    #angles_AT.append([self.Atoms[v1-1].Type,self.Atoms[key-1].Type,self.Atoms[v2-1].Type])
       
        tmp = torsions
        torsions = CheckIfUniqueListAndReturn(torsions)
        if torsions != tmp:
            sys.stderr.write("--- !Warning in "+inspect.stack()[0][3]+": torsions had double accounting. Sth went wrong. Continue.\n")
        
        tmp = angles
        angles = CheckIfUniqueListAndReturn(angles)
        if angles != tmp:
            sys.stderr.write("--- !Warning in "+inspect.stack()[0][3]+": angles had double accounting. Sth went wrong. Continue.\n")

        ### Torsions
        if self.Torsions != []:
            self.Torsions = []
            sys.stderr.write("--- !Warning in "+inspect.stack()[0][3]+": Mol.Torsions are not empty! Erase them and fill new ones according to bond connectivity. Continue.\n")
        for i,Tn in enumerate(torsions):
            T = Torsion(i+1)
            T.AtNumbers = copy(Tn)
            T.AtNames = copy(torsions_AN[i])
            T.AtTypes = copy(torsions_AT[i])
            T.Description = 'generated by ffconv'
            self.Torsions.append(copy(T))


        ### Angles
        if self.Angles != []:
            self.Angles = []
            sys.stderr.write("--- !Warning in "+inspect.stack()[0][3]+": Mol.Angles are not empty! Erase them and fill new ones according to bond connectivity. Continue.\n")
        for i,An in enumerate(angles):
            A = Angle(i+1)
            A.AtNumbers = copy(An)
            A.AtNames = copy(angles_AN[i])
            A.AtTypes = copy(angles_AT[i])
            A.Description = 'generated by ffconv'
            self.Angles.append(copy(A))




    def GenExclusions(self,NREXCL):
        exclusions=[]
        
        # TODO checks if positive integer 
        #NREXCL = self.NREXCL
       
        if isnan(NREXCL):
            sys.stderr.write("--- !Error in "+inspect.stack()[0][3]+": NREXCL is NaN. Aborting.\n")
            sys.exit(1)
            
 
        if NREXCL.__class__.__name__ != 'int':
            sys.stderr.write("--- !Error in "+inspect.stack()[0][3]+": NREXCL must be positive integer. Now ["+str(NREXCL)+"]. Aborting.\n")
            sys.exit(1)
            
            
        
        connect_map = self.GenConnectMap()

        n = 1
        Dexcl = {}
        for i in connect_map.keys():
            Dexcl[i] = []
            lst = [ i ]
            n=1
            while n <= NREXCL:
                tmpD=[]
                for v in lst:
                    if v in connect_map.keys():
                        for v2 in connect_map[v]:
                            tmpD.append(v2)
                lst = tmpD
                n += 1
                Dexcl[i] += tmpD

        for k,v in Dexcl.iteritems():
            for v2 in v:
                if v2 != k: 
                    exclusions.append([ k, v2 ])

        exclusions = CheckIfUniqueListAndReturn(exclusions)

        ### Exclusions
        if self.Exclusions != []:
            self.Exclusions = []
            sys.stderr.write("--- !Warning in "+inspect.stack()[0][3]+": Mol.Exclusions are not empty! Erase them and fill new ones according to bond connectivity. Continue.\n")
        for i,v in enumerate(exclusions):
            F = Exclusion(i+1)
            F.AtNumbers = copy(v)
            F.Description = 'generated by ffconv'
            self.Exclusions.append(copy(F))
        





    def GenFFFromMolecule(self):

        def FindAtomTypesViaAtomNumbersAndReturn(self,F):
            Types=[]
            for val in F.AtNumbers:
                for F2 in self.Atoms:
                    if ( int(val) == int(F2.Number) ):
                        Types.append( copy(F2.Type) )
            
            # Checks
            if len(Types) != len(F.AtNumbers):
                sys.stderr.write("--- !Error in "+inspect.stack()[0][3]+": could not find all atomtypes based on atom numbers: AtTypes ["+str(Types)+"], AtNumbers: ["+str(F.AtNumbers)+"],  for entry ["+str(F.__dict__)+"] of class ["+str(F.__class__.__name__)+"]. Aborting.\n")
                sys.exit(1)
                        
            return Types

        def FillAtomTypesNamesBasedOnAtomNumbers(self):
            lists=['Bonds','Pairs','Angles','Torsions','Impropers']
            for l in lists:
                for F in getattr(self,l):
                    if F.AtTypes!=['none','none'] and F.AtTypes!=['none','none','none'] and F.AtTypes!=['none','none','none','none']:
                        sys.stderr.write("--- !Warning in "+inspect.stack()[0][3]+": overwriting not none entries in AtTypes ["+str(F.AtTypes)+"] for list ["+l+"]\n")
                    setattr( F,'AtTypes',FindAtomTypesViaAtomNumbersAndReturn(self,F) )
     
            
        FillAtomTypesNamesBasedOnAtomNumbers(self)
        
        # Build FF based on Mol
        FF = ForceField( self.FFName )
        FF.Units = copy(self.Units)
        FFFamily = GetFFFamily( self.FFName ) 
        FF.NREXCL = copy(FFvars[FFFamily]['NREXCL'])
        FF.GeneratedBy = inspect.stack()[0][3]
        # Mind! Some FF parameters will be filled out based on the FFName

        for F in self.Atoms:
            AT=AtomType()
            cpAllAttributes(F,AT,['Charge','X','Y','Z','Number','Name'])
            FF.AtomTypes.append( AT )
            
        for F in self.Bonds:
            T=BondType()
            cpAllAttributes(F,T,['Number','AtNames','AtNumbers'])
            FF.BondTypes.append( T )
        
        for F in self.Pairs:
            T=PairType()
            cpAllAttributes(F,T,['Number','AtNames','AtNumbers'])
            FF.PairTypes.append( T )

        for F in self.Angles:
            T=AngleType()
            cpAllAttributes(F,T,['Number','AtNames','AtNumbers'])
            FF.AngleTypes.append( T )

        for F in self.Torsions:
            T=TorsionType()
            cpAllAttributes(F,T,['Number','AtNames','AtNumbers'])
            FF.TorsionTypes.append( T )

        for F in self.Impropers:
            T=ImproperType()
            cpAllAttributes(F,T,['Number','AtNames','AtNumbers'])
            FF.ImproperTypes.append( T )
        
        # TODO: check if yopu need this. CHARMM allows identical entries for FF ?
        FF.RemoveIdenticalFFEntries() 
        return FF


    def LoadCoordinatesFromPdb(self,fn):
        
        ifpdb = open(fn,'r')        

        # READING PDB COORDINATES FROM THE PDB FILE
        M = Molecule('tmp')
        cnt=0
        for line in ifpdb:
            ln=line.split()
            if len(ln) < 3: continue
            if ln[0]=='ATOM' or ln[0]=='HETATM': 
                cnt+=1
                F = Atom(cnt)
                F.Name=ln[2]
                line2=line[31:55]
                ln2=line2.split()
                F.X=float(ln2[0])
                F.Y=float(ln2[1])
                F.Z=float(ln2[2])
                M.Atoms.append(F)                

        if len(M.Atoms) != len(self.Atoms):
            sys.stderr.write("--- !Error in "+inspect.stack()[0][3]+": number of atoms in PDB file ["+str(M.Atoms)+"] does not match the number of atoms in the molecule instance ["+str(len(self.Atoms))+"]. Aborting.\n")
            sys.exit(1)
            
        for F in self.Atoms:
            F2 = FindObjByFieldValue(M.Atoms,'Number',F.Number)
            if F2.Name != F.Name:
                sys.stderr.write("--- !Error in "+inspect.stack()[0][3]+": atom names for atoms with the same number ["+str(F.Number)+"] are different in the pdb file ["+F2.Name+"] and in the molecule instance ["+F.Name+"]. Aborting.\n")
                sys.exit(1)
            F.X = F2.X
            F.Y = F2.Y
            F.Z = F2.Z
        return


    def DihedralsFourier2RB(self):
        for F in self.Torsions + self.Impropers:
            if F.Func == 'Fourier':
                F.p = copy(Fourier2RB(F.p))
                F.Func = 'RB'
                

        

###--------------------------------------------------------------------------------------------------------------------------------
### Some auxiliary functions and general functions 

def CheckIfAllElementsInArrayAreTheSameAndReturn(Arr):
        for i,v1 in enumerate(Arr):
            for j,v2 in enumerate(Arr[i+1:]):
                if v1 != v2:
                    sys.stderr.write('--- !Warning in function '+inspect.stack()[0][3]+': not all elements in in the given array identical: e.g. ['+v1+'] and ['+v2+']. Return empty list.\n')
                    return []
        return Arr[0]


def CheckIfUniqueListAndReturn(lst):
        unique_lst = []
        for F in lst:
            if F not in unique_lst and F[::-1] not in unique_lst:
                unique_lst.append(copy(F))
        return unique_lst


def Fourier2RB(p):
    if len(p) < 6: p = p + (6-len(p))*[0.0]
    C0 = p[1] + 0.5 * (p[0] + p[2])
    C1 = 0.5 * ( - p[0] + 3.0*p[2] )
    C2 = - p[1] + 4.0 * p[3]
    C3 = - 2.0 * p[2]
    C4 = - 4.0 * p[3]
    C5 = 0.0
    pRB = [C0,C1,C2,C3,C4,C5]
    return pRB


